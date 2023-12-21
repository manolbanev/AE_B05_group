import matplotlib.pyplot as plt
from WP5.stress_states import (get_safety_margin, get_sigma_absolute, get_stringer_buckling,
                               get_skin_buckling_function, get_shear_stress, get_web_buckling, get_sigma)
from WP5.components import wing1


def plot_safety_margin():
    fig, ax = plt.subplots(3)
    x = wing1.span
    ax[0].plot(x, get_safety_margin()[0](x), color='magenta')
    ax[0].set_title('Stringer Safety Margin')
    ax[1].plot(x, get_safety_margin()[1](x), color='cyan')
    ax[1].set_title('Skin Safety Margin')
    ax[2].plot(x, get_safety_margin()[2](x), color='purple')
    ax[2].set_title('Web Safety Margin')
    plt.tight_layout()
    ax[0].set_ylim(0, 10)
    ax[1].set_ylim(0, 10)
    ax[2].set_ylim(0, 10)
    ax[0].axhline(y=1, color='black', linestyle='dashed')
    ax[1].axhline(y=1, color='black', linestyle='dashed')
    ax[2].axhline(y=1, color='black', linestyle='dashed')
    plt.show()


def plot_failures():
    figs, axs = plt.subplots(5)
    plt.tight_layout()
    axs[0].plot(wing1.span, get_sigma(), color='magenta')
    axs[0].set_title('Normal stress')
    axs[1].plot(wing1.span, get_stringer_buckling(), color='cyan')
    axs[1].set_title('Stringer Buckling')
    axs[2].plot(wing1.span, get_web_buckling(), color='lime')
    axs[2].set_title('Web Buckling')
    axs[3].plot(wing1.span, get_shear_stress(), color='yellow')
    axs[3].set_title('Shear Stress')
    axs[4].plot(wing1.span, get_skin_buckling_function(), color='purple')
    axs[4].set_title('Skin Buckling')
    plt.show()


def plot_failures_combined():
    fig, ax = plt.subplots(2)
    x = wing1.span
    ax[0].plot(x, get_sigma_absolute(), color='magenta')
    ax[0].plot(x, get_stringer_buckling(), color='cyan', linestyle='dashed')
    ax[0].plot(x, get_skin_buckling_function(), color='purple', linestyle='dashed')
    ax[0].legend(['Normal stress', 'Stringer Buckling', 'Skin Buckling'])
    ax[0].set_xlabel('Span [m]')
    ax[0].set_ylabel('Normal stress [MPa]')
    ax[1].plot(x, get_shear_stress(), color='yellow')
    ax[1].plot(x, get_web_buckling(), color='lime', linestyle='dashed')
    ax[1].legend(['Shear stress', 'Web Buckling'])
    ax[1].set_xlabel('Span [m]')
    ax[1].set_ylabel('Shear stress [MPa]')
    plt.show()

