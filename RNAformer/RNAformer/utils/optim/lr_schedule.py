import math

import numpy as np


def get_learning_rate_schedule(scheduler_config):
    def lr_lambda(current_step: int):

        training_steps = scheduler_config.num_training_steps - scheduler_config.num_warmup_steps

        if current_step < scheduler_config.num_warmup_steps:
            return float(current_step) / float(max(1, scheduler_config.num_warmup_steps))
        elif scheduler_config.schedule == 'linear':
            linear_decay = max(0.0, float(scheduler_config.num_training_steps - current_step) / float(
                max(1, training_steps)))
            return scheduler_config.decay_factor + (1 - scheduler_config.decay_factor) * linear_decay
        elif scheduler_config.schedule == 'cosine':
            cosine_decay = max(0.0, (1 + math.cos(
                math.pi * (current_step - scheduler_config.num_warmup_steps) / float(max(1, training_steps)))) / 2)
            return scheduler_config.decay_factor + (1 - scheduler_config.decay_factor) * cosine_decay
        elif scheduler_config.schedule == 'logstep':
            interval_length = training_steps // scheduler_config.num_steps
            log_space_lrs = np.logspace(np.log10(1.0), np.log10(1.0 * scheduler_config.decay_factor),
                                        scheduler_config.num_steps)
            current_interval = (current_step - scheduler_config.num_warmup_steps) // interval_length
            current_interval = min(current_interval, scheduler_config.num_steps - 1)
            return log_space_lrs[current_interval]
        elif scheduler_config.schedule == 'step':
            interval_length = training_steps // scheduler_config.num_steps
            lr_decrement = (1.0 - 1.0 * scheduler_config.decay_factor) / (scheduler_config.num_steps - 1)
            linear_space_lrs = [1.0 - lr_decrement * i for i in range(scheduler_config.num_steps)]
            current_interval = (current_step - scheduler_config.num_warmup_steps) // interval_length
            current_interval = min(current_interval, scheduler_config.num_steps - 1)
            return linear_space_lrs[current_interval]
        elif scheduler_config.schedule == 'const':
            return 1.0
        elif scheduler_config.schedule == 'wsd':
            if current_step > (scheduler_config.num_training_steps - scheduler_config.num_warmup_steps):
                return float(scheduler_config.num_training_steps - current_step) / float(
                    max(1, scheduler_config.num_warmup_steps))
            else:
                return 1.0

    return lr_lambda


if __name__ == "__main__":
    from argparse import Namespace
    import matplotlib.pyplot as plt

    num_training_steps = 50000

    scheduler_step = Namespace(num_warmup_steps=1000, num_training_steps=num_training_steps, decay_factor=0.001,
                               num_steps=5, schedule='step')
    scheduler_logstep = Namespace(num_warmup_steps=1000, num_training_steps=num_training_steps, decay_factor=0.001,
                                  num_steps=5, schedule='logstep')
    scheduler_linear = Namespace(num_warmup_steps=1000, num_training_steps=num_training_steps, decay_factor=0.001,
                                 num_steps=5, schedule='linear')
    scheduler_cosine = Namespace(num_warmup_steps=1000, num_training_steps=num_training_steps, decay_factor=0.001,
                                 num_steps=5, schedule='cosine')
    scheduler_const = Namespace(num_warmup_steps=1000, num_training_steps=num_training_steps, decay_factor=0.001,
                                num_steps=5, schedule='const')

    lr_step = get_learning_rate_schedule(scheduler_step)
    lr_logstep = get_learning_rate_schedule(scheduler_logstep)
    lr_linear = get_learning_rate_schedule(scheduler_linear)
    lr_cosine = get_learning_rate_schedule(scheduler_cosine)
    lr_const = get_learning_rate_schedule(scheduler_const)

    step_schedule = [lr_step(i) for i in range(num_training_steps)]
    logstep_schedule = [lr_logstep(i) for i in range(num_training_steps)]
    linear_schedule = [lr_linear(i) for i in range(num_training_steps)]
    cosine_schedule = [lr_cosine(i) for i in range(num_training_steps)]
    const_schedule = [lr_const(i) for i in range(num_training_steps)]

    plt.plot(step_schedule, label='step')
    plt.plot(logstep_schedule, label='logstep')
    plt.plot(linear_schedule, label='linear')
    plt.plot(cosine_schedule, label='cosine')
    plt.plot(const_schedule, label='const')
    plt.legend()
    plt.show()
