# sitecustomize.py
"""
Force a PyTorch program to run entirely on CPU without editing the codebase.

- Hides CUDA/MPS from PyTorch (cuda.is_available() -> False, device_count -> 0)
- Redirects .cuda() and .to('cuda') calls on Tensors and Modules to CPU
- Maps torch.load / torch.jit.load checkpoints to CPU by default
- Forces safetensors loads to CPU if available
- Removes 'cache_key' kwarg if forward() doesn't accept it (API mismatch guard)
- Provides CPU-safe stubs for flash-attn:
    * flash_attn.bert_padding: unpad_input / pad_input
    * flash_attn.flash_attn_interface: flash_attn_varlen_qkvpacked_func
  and mirrors them to builtins for bare-name calls.
"""

# =============================================================================
# flash-attn CPU FALLBACKS (define BEFORE importing torch-dependent modules)
# =============================================================================
# -- bert_padding: unpad_input / pad_input ------------------------------------
try:
    import importlib, types, sys, builtins
    # ensure package skeleton exists so later imports don't fail
    if "flash_attn" not in sys.modules:
        pkg = types.ModuleType("flash_attn"); pkg.__path__ = []
        sys.modules["flash_attn"] = pkg

    def _normalize_keep(mask):
        import torch
        if mask is None:
            return None
        if mask.dtype == torch.bool:
            # Heuristic: True means PAD in many libs -> keep = ~mask
            return ~mask if mask.any() else mask
        return mask != 0

    def _cpu_unpad_input(x, key_padding_mask=None):
        import torch
        assert x.dim() == 3, f"expected [B,L,D], got {tuple(x.shape)}"
        B, L, D = x.shape
        device = x.device
        keep = _normalize_keep(key_padding_mask)

        if keep is None:
            x_unpad = x.reshape(B * L, D)
            indices = torch.arange(B * L, device=device, dtype=torch.long)
            seqlens = torch.full((B,), L, device=device, dtype=torch.int32)
            cu_seqlens = torch.zeros(B + 1, device=device, dtype=torch.int32)
            cu_seqlens[1:] = seqlens.cumsum(0)
            max_s = int(L)
            return x_unpad, indices, cu_seqlens, max_s

        keep = keep.to(device=device, dtype=torch.bool)
        seqlens = keep.sum(dim=1).to(torch.int32)
        cu_seqlens = torch.zeros(B + 1, device=device, dtype=torch.int32)
        cu_seqlens[1:] = seqlens.cumsum(0)
        max_s = int(seqlens.max().item()) if seqlens.numel() else 0

        x_flat = x.reshape(B * L, D)
        lin = torch.arange(B * L, device=device, dtype=torch.long).reshape(B, L)
        indices = lin[keep]                      # [nnz]
        x_unpad = x_flat[indices]                # [nnz, D]
        return x_unpad, indices, cu_seqlens, max_s

    def _as_int_maybe(x):
        """Try to coerce scalar-like input (int, 0-dim tensor/np) to int; else return None."""
        try:
            # torch / numpy scalar
            if hasattr(x, "ndim") and int(getattr(x, "ndim")) == 0:
                return int(x.item())
        except Exception:
            pass
        try:
            # plain python int-ish
            return int(x)
        except Exception:
            return None

    def _cpu_pad_input(x_unpad, indices=None, batch=None, seqlens=None, L=None):
        """
        CPU fallback for flash-attn pad_input.
        Supports:
          - pad_input(x_unpad, indices, L=L)
          - pad_input(x_unpad, indices, batch=B, seqlens=seqlens)
            where seqlens is:
              * scalar (max length L), or
              * per-sequence lengths: shape [B], or
              * cumulative seqlens (cu_seqlens): shape [B+1]
        Returns [B, L, D].
        """
        import torch
        assert indices is not None, "pad_input needs 'indices'"

        # If seqlens is a scalar, treat it as L
        if L is None and seqlens is not None:
            L_scalar = _as_int_maybe(seqlens)
            if L_scalar is not None:
                L = L_scalar
            else:
                seqlens = torch.as_tensor(seqlens, device=indices.device)
                if seqlens.dim() == 1:
                    if seqlens.numel() >= 2 and torch.all(seqlens[1:] >= seqlens[:-1]) and seqlens[0].item() == 0:
                        # cu_seqlens
                        B = int(seqlens.numel() - 1)
                        lens = (seqlens[1:] - seqlens[:-1])
                        L = int(lens.max().item()) if lens.numel() else 0
                    else:
                        # per-seq lengths
                        B = int(seqlens.numel())
                        L = int(seqlens.max().item()) if seqlens.numel() else 0
                else:
                    # Non-1D non-scalar -> fallback
                    raise AssertionError("seqlens must be 1D or scalar when provided")

        if L is None:
            assert batch is not None, "pad_input needs either L or (batch & seqlens)"
            B = int(batch)
            max_index = int(indices.max().item()) if indices.numel() else -1
            L = (max_index + 1 + (B - 1)) // B if B > 0 else 0
        else:
            if batch is None:
                max_index = int(indices.max().item()) if indices.numel() else -1
                B = (max_index + 1) // int(L) + 1 if L else 1
            else:
                B = int(batch)

        L = int(L)
        D = x_unpad.shape[-1]
        device = x_unpad.device

        out = torch.zeros(B * L, D, device=device, dtype=x_unpad.dtype)
        if indices.numel():
            out[indices] = x_unpad
        return out.reshape(B, L, D)

    # Install submodule with stubs if real module is missing
    try:
        importlib.import_module("flash_attn.bert_padding")
    except Exception:
        mod_bp = types.ModuleType("flash_attn.bert_padding")
        mod_bp.unpad_input = _cpu_unpad_input
        mod_bp.pad_input = _cpu_pad_input
        sys.modules["flash_attn.bert_padding"] = mod_bp

    # ALWAYS mirror to builtins so bare-name calls work
    builtins.unpad_input = _cpu_unpad_input   # type: ignore[attr-defined]
    builtins.pad_input   = _cpu_pad_input     # type: ignore[attr-defined]

except Exception:
    pass

# -- flash_attn_interface: flash_attn_varlen_qkvpacked_func -------------------
try:
    import importlib, types, sys, builtins
    import torch
    import torch.nn.functional as F

    def _flash_stub_varlen_qkvpacked(qkv, cu_seqlens, max_s: int,
                                     dropout_p: float = 0.0,
                                     softmax_scale = None,
                                     causal: bool = False,
                                     window_size = None,
                                     return_softmax: bool = False,
                                     **kwargs):
        """
        CPU fallback approximating flash attention for varlen-packed QKV.
        qkv: [nnz, 3, h, d]; cu_seqlens: [B+1] (int32); returns out: [nnz, h, d]
        """
        assert qkv.dim() == 4 and qkv.size(1) == 3, f"expected [nnz,3,h,d], got {tuple(qkv.shape)}"
        nnz, _, h, d = qkv.shape
        device = qkv.device
        dtype  = qkv.dtype

        out = torch.empty((nnz, h, d), device=device, dtype=dtype)
        B = int(cu_seqlens.numel() - 1)
        scale = (d ** -0.5) if softmax_scale is None else float(softmax_scale)

        start = 0
        for b in range(B):
            s = int((cu_seqlens[b+1] - cu_seqlens[b]).item())
            if s <= 0:
                continue
            # slice to [s, 3, h, d] -> [3, h, s, d]
            q, k, v = qkv[start:start+s].permute(1, 2, 0, 3).unbind(0)  # each: [h, s, d]
            attn = (q @ k.transpose(-1, -2)) * scale                    # [h, s, s]
            if causal:
                causal_mask = torch.ones((s, s), device=device, dtype=torch.bool).triu(1)
                attn.masked_fill_(causal_mask, float("-inf"))
            attn = F.softmax(attn, dim=-1)
            out_b = attn @ v                                            # [h, s, d]
            out[start:start+s] = out_b.permute(1, 0, 2)                 # -> [s, h, d]
            start += s

        if return_softmax:
            return out, None
        return out

    try:
        importlib.import_module("flash_attn.flash_attn_interface")
    except Exception:
        if "flash_attn" not in sys.modules:
            pkg = types.ModuleType("flash_attn"); pkg.__path__ = []
            sys.modules["flash_attn"] = pkg
        mod_fi = types.ModuleType("flash_attn.flash_attn_interface")
        mod_fi.flash_attn_varlen_qkvpacked_func = _flash_stub_varlen_qkvpacked
        sys.modules["flash_attn.flash_attn_interface"] = mod_fi
    else:
        sys.modules["flash_attn.flash_attn_interface"].flash_attn_varlen_qkvpacked_func = _flash_stub_varlen_qkvpacked

    builtins.flash_attn_varlen_qkvpacked_func = _flash_stub_varlen_qkvpacked  # type: ignore[attr-defined]
except Exception:
    pass

# =============================================================================
# TORCH CPU FORCING & SHIMS
# =============================================================================
try:
    import torch
    from typing import Any
    import inspect

    # Hide accelerators
    try:
        torch.cuda.is_available = lambda: False  # type: ignore[attr-defined,assignment]
        torch.cuda.device_count = lambda: 0      # type: ignore[attr-defined,assignment]
    except Exception:
        pass

    # Prefer CPU default device (PyTorch ≥ 2.0)
    try:
        torch.set_default_device("cpu")  # type: ignore[attr-defined]
    except Exception:
        pass

    # Helpers
    def _force_cpu_device(dev: Any) -> Any:
        try:
            if dev is None:
                return dev
            s = str(dev).lower()
            if s.startswith("cuda") or s.startswith("mps"):
                return torch.device("cpu")
        except Exception:
            pass
        return dev

    # Redirect .cuda() to .cpu()
    try:
        def _to_cpu(self, *args, **kwargs):
            return self.cpu()
        torch.Tensor.cuda = _to_cpu           # type: ignore[assignment]
        torch.nn.Module.cuda = _to_cpu        # type: ignore[assignment]
    except Exception:
        pass

    # Coerce .to("cuda") to CPU (Tensor)
    try:
        _orig_tensor_to = torch.Tensor.to
        def _tensor_to(self, *args, **kwargs):
            if args:
                a = list(args)
                a[0] = _force_cpu_device(a[0])
                args = tuple(a)
            if "device" in kwargs:
                kwargs["device"] = _force_cpu_device(kwargs["device"])
            return _orig_tensor_to(self, *args, **kwargs)
        torch.Tensor.to = _tensor_to  # type: ignore[assignment]
    except Exception:
        pass

    # Coerce .to("cuda") to CPU (Module)
    try:
        _orig_module_to = torch.nn.Module.to
        def _module_to(self, *args, **kwargs):
            if args:
                a = list(args)
                a[0] = _force_cpu_device(a[0])
                args = tuple(a)
            if "device" in kwargs:
                kwargs["device"] = _force_cpu_device(kwargs["device"])
            return _orig_module_to(self, *args, **kwargs)
        torch.nn.Module.to = _module_to  # type: ignore[assignment]
    except Exception:
        pass

    # Force checkpoint loads to CPU
    try:
        _orig_load = torch.load
        def _load_cpu(*args, **kwargs):
            if "map_location" not in kwargs or kwargs["map_location"] is None:
                kwargs["map_location"] = "cpu"
            return _orig_load(*args, **kwargs)
        torch.load = _load_cpu  # type: ignore[assignment]
    except Exception:
        pass

    try:
        _orig_jit_load = torch.jit.load  # type: ignore[attr-defined]
        def _jit_load_cpu(*args, **kwargs):
            if "map_location" not in kwargs or kwargs["map_location"] is None:
                kwargs["map_location"] = "cpu"
            return _orig_jit_load(*args, **kwargs)
        torch.jit.load = _jit_load_cpu  # type: ignore[assignment]
    except Exception:
        pass

    # Safetensors to CPU
    try:
        from safetensors.torch import load_file as _sf_load_file  # type: ignore
        def _sf_load_file_cpu(*args, **kwargs):
            kwargs["device"] = "cpu"
            return _sf_load_file(*args, **kwargs)
        import safetensors.torch as _sft  # type: ignore
        _sft.load_file = _sf_load_file_cpu  # type: ignore[assignment]
    except Exception:
        pass

    # Strip 'cache_key' if forward() doesn't accept it
    try:
        _orig_call_impl = torch.nn.Module._call_impl
        def _call_impl_strip_cache_key(self, *args, **kwargs):
            if "cache_key" in kwargs:
                try:
                    sig = inspect.signature(self.forward)
                    if "cache_key" not in sig.parameters:
                        kwargs = dict(kwargs)
                        kwargs.pop("cache_key", None)
                except Exception:
                    pass
            return _orig_call_impl(self, *args, **kwargs)
        torch.nn.Module._call_impl = _call_impl_strip_cache_key  # type: ignore[assignment]
    except Exception:
        pass

except Exception:
    # Never break startup
    pass
