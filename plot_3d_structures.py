#!/usr/bin/env python3
import os
import json
import time
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import urllib.request
from Bio.PDB import MMCIFParser, PDBParser  # future use


# --------------------------- Utilities ---------------------------

def download_pdb(pdb_id, save_dir="evaluation/predictions/gt", file_format="cif"):
    os.makedirs(save_dir, exist_ok=True)
    dst = Path(save_dir) / f"{pdb_id}_gt.{file_format}"
    if dst.exists():
        print(f"⚡ GT structure already exists: {dst}")
        return dst
    url = f"https://files.rcsb.org/download/{pdb_id}.{file_format}"
    print(f"⬇️ Downloading GT structure {pdb_id} to {dst}")
    urllib.request.urlretrieve(url, dst)
    print(f"✅ Download complete: {dst}")
    return dst


def _run_pymol(script: str):
    td = tempfile.mkdtemp()
    pml = Path(td) / "render.pml"
    pml.write_text(script, encoding='utf8')
    print(f"🔧 Running PyMOL script: {pml}")

    proc = subprocess.run(["pymol", "-cq", str(pml)], capture_output=True, text=True)
    if proc.returncode != 0:
        print(f">>> PyMOL returned code {proc.returncode}")
        print(proc.stdout)
        print(proc.stderr)
    else:
        print(f"✅ PyMOL finished successfully for script {pml}")
    try:
        pml.unlink()
        shutil.rmtree(td)
    except Exception:
        pass


def _python_block_apply_canonical_view(pdb_id: str, views_dir: Path, which: str,
                                       buffer: float, transparent: bool,
                                       force_rebuild: bool):
    """
    Return a PML 'python' block that:
      - if a cached view exists, loads and applies it;
      - otherwise (or if force), computes from GT-only selection,
        saves it, and applies it.

    which: "complex" uses selection 'gt'
           "rna"     uses selection 'gt_rna'
    """
    views_dir = Path(views_dir)
    views_dir.mkdir(parents=True, exist_ok=True)
    view_file = views_dir / f"{pdb_id.upper()}_{which}.json"

    # This block runs inside PyMOL's Python.
    return f"""
python
from pymol import cmd
import json, os
pdb_id = {json.dumps(pdb_id.upper())}
which = {json.dumps(which)}
view_path = {json.dumps(str(view_file))}
sel = 'gt_rna' if which == 'rna' else 'gt'
# Always base camera on GT-only selection
cmd.set('auto_zoom', 0)

rebuild = {str(bool(force_rebuild))}
if (not rebuild) and os.path.exists(view_path):
    # Load and apply cached canonical view
    try:
        with open(view_path, 'r') as f:
            v = tuple(json.load(f))
        cmd.set_view(v)
        print(f"📥 Applied cached {{which}} canonical view for {{pdb_id}} from: {{view_path}}")
    except Exception as e:
        print(f"⚠️ Failed to read cached view {{view_path}}: {{e}}; will rebuild.")
        rebuild = True

if rebuild or (not os.path.exists(view_path)):
    # Compute canonical view from GT selection only
    cmd.orient(sel)
    cmd.zoom(sel, buffer={float(buffer)}, complete=1)
    cmd.clip('all', 0)
    v = tuple(cmd.get_view())
    # Save and apply
    try:
        with open(view_path, 'w') as f:
            json.dump(list(v), f)
        print(f"💾 Saved canonical {{which}} view for {{pdb_id}} to: {{view_path}}")
    except Exception as e:
        print(f"⚠️ Could not save view to {{view_path}}: {{e}}")
    cmd.set_view(v)
python end
"""


# --------------------------- Renderers ---------------------------

def render_complex(
    file_path,
    gt_file,
    output_image,
    pred_color,
    gt_color,
    pred_rna_color='yellow',
    gt_rna_color='palecyan',
    *,
    buffer=20,
    orthographic=True,
    width=1000, height=1000,
    dpi=300,
    ray=1,
    shadows=False,
    transparent=False,
    views_dir="pymol_views",
    pdb_id="XXXX",
    force_view_cache=False
):
    """
    Render full complex view. Camera is locked to GT-only canonical view
    (per PDB, cached in views_dir).
    """
    output_image = Path(output_image).absolute()
    if output_image.exists():
        print(f"⚡ Skipping complex render; file exists: {output_image}")
        return
    output_image.parent.mkdir(exist_ok=True, parents=True)

    bg = "transparent" if transparent else "white"
    ray_bg = "off" if transparent else "on"

    script = f'''
reinitialize
set auto_zoom, off
set orthoscopic, {"on" if orthographic else "off"}
set depth_cue, off
set ray_shadows, {"on" if shadows else "off"}
set ray_opaque_background, {ray_bg}
set max_ups, 0
set antialias, 2

load {gt_file}, gt
load {file_path}, pred

# Align prediction to GT (on all atoms) without changing camera
align pred, gt

# Build RNA selections (used for coloring; camera is GT based)
select pred_rna, pred and polymer.nucleic
select gt_rna,   gt   and polymer.nucleic

# Apply canonical camera based on GT only (cached)
{_python_block_apply_canonical_view(pdb_id, Path(views_dir), "complex", buffer, transparent, force_view_cache)}

# Show complex
hide everything
show cartoon, pred
color {pred_color}, pred
show lines, gt
color {gt_color}, gt

# Emphasize RNA parts
show cartoon, pred_rna
color {pred_rna_color}, pred_rna
show cartoon, gt_rna
color {gt_rna_color}, gt_rna

bg_color {bg}
png {output_image}, width={width}, height={height}, dpi={dpi}, ray={ray}
quit
'''
    _run_pymol(script)
    if output_image.exists():
        print(f"✅ Complex image created: {output_image}")
    else:
        print(f"❌ Failed to create complex image: {output_image}")


def render_rna_only(
    file_path,
    gt_file,
    output_image,
    pred_rna_color,
    gt_rna_color,
    *,
    buffer=20,
    orthographic=True,
    width=1000, height=1000,
    dpi=300,
    ray=1,
    shadows=False,
    transparent=False,
    views_dir="pymol_views",
    pdb_id="XXXX",
    force_view_cache=False
):
    """
    Render RNA-only view. Camera locked to GT-RNA canonical view (per PDB).
    """
    output_image = Path(output_image).absolute()
    if output_image.exists():
        print(f"⚡ Skipping RNA-only render; file exists: {output_image}")
        return
    output_image.parent.mkdir(exist_ok=True, parents=True)

    bg = "transparent" if transparent else "white"
    ray_bg = "off" if transparent else "on"

    script = f'''
reinitialize
set auto_zoom, off
set orthoscopic, {"on" if orthographic else "off"}
set depth_cue, off
set ray_shadows, {"on" if shadows else "off"}
set ray_opaque_background, {ray_bg}
set max_ups, 0
set antialias, 2

load {gt_file}, gt
load {file_path}, pred

# RNA selections and alignment on RNA only
select gt_rna,   gt   and polymer.nucleic
select pred_rna, pred and polymer.nucleic
align pred_rna, gt_rna

# Apply canonical camera based on GT RNA only (cached)
{_python_block_apply_canonical_view(pdb_id, Path(views_dir), "rna", buffer, transparent, force_view_cache)}

# Show only RNA
hide everything
show cartoon, pred_rna
color {pred_rna_color}, pred_rna
show cartoon, gt_rna
color {gt_rna_color}, gt_rna

bg_color {bg}
png {output_image}, width={width}, height={height}, dpi={dpi}, ray={ray}
quit
'''
    _run_pymol(script)
    if output_image.exists():
        print(f"✅ RNA-only image created: {output_image}")
    else:
        print(f"❌ Failed to create RNA-only image: {output_image}")


# --------------------------- Combining images ---------------------------

def combine_side_by_side(img1, img2, out_png, title1, title2, show):
    out_png = Path(out_png)
    if out_png.exists():
        print(f"⚡ Skipping combine; file exists: {out_png}")
        return
    if not (Path(img1).exists() and Path(img2).exists()):
        missing = [f for f in (img1, img2) if not Path(f).exists()]
        print(f"❌ Cannot combine; missing images: {missing}")
        return
    print(f"🔗 Combining images:\n  1: {img1}\n  2: {img2}\n  into: {out_png}")
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax, img, title in zip(axes, [img1, img2], [title1, title2]):
        ax.imshow(mpimg.imread(img))
        ax.set_title(title)
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print(f"✅ Combined image saved: {out_png}")
    if show:
        plt.show()
    plt.close(fig)


# --------------------------- Driver ---------------------------

def main(
    rnaformer_pred_dir,
    alphafold_pred_dir,
    gt_dir,
    out_dir,
    show_plot=False,
    specific_id=None,
    render_rf_only=False,
    plot_rna_only=False,
    rf_color='red',
    af3_color='blue',
    gt_color='gray',
    gt_rna_color='palecyan',
    pred_rna_color='yellow',
    alg1_name='AF3',
    alg2_name='DSSR',
    *,
    buffer=20,
    orthographic=True,
    img_width=1000,
    img_height=1000,
    dpi=300,
    ray=1,
    shadows=False,
    transparent=False,
    views_dir="pymol_views",
    force_view_cache=False
):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    views_dir = Path(views_dir)
    views_dir.mkdir(exist_ok=True, parents=True)

    print(
        "🏁 Starting processing:\n"
        f"  rnaformer_pred_dir = {rnaformer_pred_dir}\n"
        f"  alphafold_pred_dir = {alphafold_pred_dir}\n"
        f"  gt_dir = {gt_dir}\n"
        f"  out_dir = {out_dir}\n"
        f"  views_dir = {views_dir}\n"
        f"  specific_id = {specific_id}\n"
        f"  buffer = {buffer}, orthographic = {orthographic}, size = {img_width}x{img_height}, dpi = {dpi}, ray = {ray}, shadows = {shadows}, transparent = {transparent}\n"
        f"  force_view_cache = {force_view_cache}"
    )

    for cif in Path(rnaformer_pred_dir).glob("*.cif"):
        pdb_id = cif.stem[:4].upper()
        if specific_id and pdb_id != specific_id.upper():
            continue
        print(f"\n⏳ Processing PDB ID: {pdb_id}")

        gt_cif = Path(gt_dir) / f"{pdb_id}_gt.cif"
        if not gt_cif.exists():
            gt_pdb = Path(gt_dir) / f"{pdb_id}.pdb"
            gt_cif = gt_pdb if gt_pdb.exists() else download_pdb(pdb_id, save_dir=gt_dir)

        # Try both lowercase & uppercase patterns for AF3/RNAfold dir
        af3_list = list(Path(alphafold_pred_dir).glob(f"*{pdb_id.lower()}*.cif")) + \
                   list(Path(alphafold_pred_dir).glob(f"*{pdb_id.upper()}*.cif"))
        af3 = af3_list[0] if af3_list else None
        rf = cif

        af3_cpx = out_dir / f"{pdb_id}_{alg1_name}_rnaprotein.png"
        rf_cpx  = out_dir / f"{pdb_id}_{alg2_name}_rnaprotein.png"
        af3_rna = out_dir / f"{pdb_id}_{alg1_name}_rna_rnaprotein.png"
        rf_rna  = out_dir / f"{pdb_id}_{alg2_name}_rna_rnaprotein.png"

        try:
            if not render_rf_only and af3:
                print(f"🔍 Rendering {alg1_name} complex for {pdb_id}")
                render_complex(
                    af3, gt_cif, af3_cpx, af3_color, gt_color,
                    pred_rna_color, gt_rna_color,
                    buffer=buffer, orthographic=orthographic,
                    width=img_width, height=img_height,
                    dpi=dpi, ray=ray, shadows=shadows,
                    transparent=transparent,
                    views_dir=views_dir, pdb_id=pdb_id,
                    force_view_cache=force_view_cache
                )
        except Exception as e:
            print(f"❌ Error rendering {alg1_name} complex for {pdb_id}: {e}")

        try:
            print(f"🔍 Rendering {alg2_name} complex for {pdb_id}")
            render_complex(
                rf, gt_cif, rf_cpx, rf_color, gt_color,
                pred_rna_color, gt_rna_color,
                buffer=buffer, orthographic=orthographic,
                width=img_width, height=img_height,
                dpi=dpi, ray=ray, shadows=shadows,
                transparent=transparent,
                views_dir=views_dir, pdb_id=pdb_id,
                force_view_cache=force_view_cache
            )
        except Exception as e:
            print(f"❌ Error rendering {alg2_name} complex for {pdb_id}: {e}")

        try:
            if not render_rf_only and plot_rna_only and af3:
                print(f"🔍 Rendering {alg1_name} RNA-only for {pdb_id}")
                render_rna_only(
                    af3, gt_cif, af3_rna, af3_color, gt_rna_color,
                    buffer=buffer, orthographic=orthographic,
                    width=img_width, height=img_height,
                    dpi=dpi, ray=ray, shadows=shadows,
                    transparent=transparent,
                    views_dir=views_dir, pdb_id=pdb_id,
                    force_view_cache=force_view_cache
                )
        except Exception as e:
            print(f"❌ Error rendering {alg1_name} RNA-only for {pdb_id}: {e}")

        try:
            if plot_rna_only:
                print(f"🔍 Rendering {alg2_name} RNA-only for {pdb_id}")
                render_rna_only(
                    rf, gt_cif, rf_rna, rf_color, gt_rna_color,
                    buffer=buffer, orthographic=orthographic,
                    width=img_width, height=img_height,
                    dpi=dpi, ray=ray, shadows=shadows,
                    transparent=transparent,
                    views_dir=views_dir, pdb_id=pdb_id,
                    force_view_cache=force_view_cache
                )
        except Exception as e:
            print(f"❌ Error rendering {alg2_name} RNA-only for {pdb_id}: {e}")

        # Combine (will skip if any image missing)
        combine_side_by_side(
            af3_cpx, rf_cpx,
            out_dir / f"{pdb_id}_{alg1_name.lower()}_{alg2_name.lower()}_complex_comparison.png",
            f"{alg1_name} complex", f"{alg2_name} complex",
            show_plot
        )
        combine_side_by_side(
            af3_rna, rf_rna,
            out_dir / f"{pdb_id}_{alg1_name.lower()}_{alg2_name.lower()}_rna_comparison.png",
            f"{alg1_name} RNA only", f"{alg2_name} RNA only",
            show_plot
        )

        time.sleep(0.2)  # small breather


# --------------------------- CLI ---------------------------

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Render RNA-protein complexes with canonical, GT-locked camera')
    parser.add_argument('--rnaformer-pred-dir', required=True, help='Directory with RNAformer .cif files')
    parser.add_argument('--alphafold-pred-dir', required=True, help='Directory with AlphaFold3/RNAfold .cif files')
    parser.add_argument('--gt-dir', required=True, help='Directory with ground-truth .cif/.pdb files or where to download')
    parser.add_argument('--out-dir', required=True, help='Output directory for images')
    parser.add_argument('--views-dir', default='pymol_views', help='Directory to cache per-PDB canonical views')
    parser.add_argument('--show-plot', action='store_true', help='Display combined images interactively')
    parser.add_argument('--specific-id', type=str, help='Only process a specific 4-letter PDB ID (case-insensitive)')
    parser.add_argument('--render-rf-only', action='store_true', help='Skip rendering AF3 complex/RNA views')
    parser.add_argument('--plot-rna-only', action='store_true', help='Also render only the RNA portions')
    parser.add_argument('--rf-color', default='red', help='Color for RNAformer predictions')
    parser.add_argument('--af3-color', default='blue', help='Color for AlphaFold3/RNAfold predictions')
    parser.add_argument('--gt-color', default='gray', help='Color for GT structures in complex')
    parser.add_argument('--gt-rna-color', default='palecyan', help='Color for GT RNA-only')
    parser.add_argument('--pred-rna-color', default='yellow', help='Color for predicted RNA-only')
    parser.add_argument('--alg1-name', default='AF3', help='Label for Algorithm 1 (e.g., AF3)')
    parser.add_argument('--alg2-name', default='DSSR', help='Label for Algorithm 2 (e.g., DSSR)')

    # Camera / image controls
    parser.add_argument('--buffer', type=float, default=20.0, help='Zoom buffer (Å) around GT to avoid cropping')
    parser.add_argument('--no-orthographic', dest='orthographic', action='store_false', help='Use perspective projection instead of orthographic')
    parser.set_defaults(orthographic=True)
    parser.add_argument('--img-width', type=int, default=1000, help='Output image width in pixels')
    parser.add_argument('--img-height', type=int, default=1000, help='Output image height in pixels')
    parser.add_argument('--dpi', type=int, default=300, help='Output DPI for saved PNG')
    parser.add_argument('--ray', type=int, default=1, choices=[0, 1], help='Ray tracing: 1 for quality, 0 for speed')
    parser.add_argument('--shadows', action='store_true', help='Enable ray shadows (off by default)')

    # View cache control
    parser.add_argument('--transparent', action='store_true', help='Render PNGs with transparent background')
    parser.add_argument('--force-view-cache', action='store_true',
                        help='Recompute and overwrite the cached canonical view for the PDB ID')

    args = parser.parse_args()

    main(
        rnaformer_pred_dir=args.rnaformer_pred_dir,
        alphafold_pred_dir=args.alphafold_pred_dir,
        gt_dir=args.gt_dir,
        out_dir=args.out_dir,
        show_plot=args.show_plot,
        specific_id=args.specific_id,
        render_rf_only=args.render_rf_only,
        plot_rna_only=args.plot_rna_only,
        rf_color=args.rf_color,
        af3_color=args.af3_color,
        gt_color=args.gt_color,
        gt_rna_color=args.gt_rna_color,
        pred_rna_color=args.pred_rna_color,
        alg1_name=args.alg1_name,
        alg2_name=args.alg2_name,
        buffer=args.buffer,
        orthographic=args.orthographic,
        img_width=args.img_width,
        img_height=args.img_height,
        dpi=args.dpi,
        ray=args.ray,
        shadows=args.shadows,
        transparent=args.transparent,
        views_dir=args.views_dir,
        force_view_cache=args.force_view_cache,
    )
