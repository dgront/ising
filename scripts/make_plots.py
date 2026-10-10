import glob
import re
import numpy as np
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

L = 1024
N_LAST = 1000
N_SPINS = L * L

FILE_PATTERN = f"en_{L}_*"


# ------------------------------------------------------------
# Read files
# ------------------------------------------------------------

files = glob.glob(FILE_PATTERN)

data = []

for filename in files:

    m = re.search(rf"en_{L}_([0-9.]+)$", filename)

    if m is None:
        print(f"Skipping: {filename}")
        continue

    T = float(m.group(1))

    arr = np.loadtxt(filename)

    frame = arr[:, 0]
    energy = arr[:, 1]       # energy per spin
    acceptance = arr[:, 2]   # acceptance fraction

    data.append((T, frame, energy, acceptance, filename))


data.sort(key=lambda x: x[0])

if not data:
    raise RuntimeError("No input files found")


# ------------------------------------------------------------
# Energy vs time
# ------------------------------------------------------------

plt.figure(figsize=(10, 6))

for T, frame, energy, acceptance, filename in data:
    plt.plot(frame, energy, label=f"T={T:g}", linewidth=1)

plt.xlabel("Frame")
plt.ylabel("Energy per spin")
plt.title("Energy vs time")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("energy_vs_time.png", dpi=200)
plt.show()


# ------------------------------------------------------------
# Acceptance fraction vs time
# ------------------------------------------------------------

plt.figure(figsize=(10, 6))

for T, frame, energy, acceptance, filename in data:
    plt.plot(frame, acceptance, label=f"T={T:g}", linewidth=1)

plt.xlabel("Frame")
plt.ylabel("Acceptance fraction")
plt.title("Metropolis acceptance fraction vs time")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("acceptance_vs_time.png", dpi=200)
plt.show()


# ------------------------------------------------------------
# Statistics from last N frames
# ------------------------------------------------------------

temperatures = []
mean_E = []
var_E = []
Cv = []
mean_acceptance = []

for T, frame, energy, acceptance, filename in data:

    if len(energy) < N_LAST:
        print(
            f"Warning: {filename} contains only {len(energy)} points; "
            f"using all of them."
        )
        e_sample = energy
        a_sample = acceptance
    else:
        e_sample = energy[-N_LAST:]
        a_sample = acceptance[-N_LAST:]

    e_mean = np.mean(e_sample)

    # variance of energy per spin
    e_var = np.mean(e_sample**2) - e_mean**2

    # heat capacity per spin:
    # Cv/N = N * Var(e) / T^2
    cv = N_SPINS * e_var / T**2

    a_mean = np.mean(a_sample)

    temperatures.append(T)
    mean_E.append(e_mean)
    var_E.append(e_var)
    Cv.append(cv)
    mean_acceptance.append(a_mean)

    print(
        f"T = {T:7.4f}   "
        f"<e> = {e_mean:10.6f}   "
        f"Var(e) = {e_var:12.6e}   "
        f"Cv/N = {cv:10.6f}   "
        f"<acc> = {a_mean:8.4f}"
    )


temperatures = np.asarray(temperatures)
mean_E = np.asarray(mean_E)
var_E = np.asarray(var_E)
Cv = np.asarray(Cv)
mean_acceptance = np.asarray(mean_acceptance)


# ------------------------------------------------------------
# E(T)
# ------------------------------------------------------------

plt.figure(figsize=(7, 5))

plt.plot(
    temperatures,
    mean_E,
    "o-"
)

plt.xlabel("Temperature T")
plt.ylabel(r"$\langle e \rangle$")
plt.title(f"Mean energy per spin, last {N_LAST} frames")
plt.tight_layout()
plt.savefig("energy_vs_temperature.png", dpi=200)
plt.show()


# ------------------------------------------------------------
# Cv(T)
# ------------------------------------------------------------

plt.figure(figsize=(7, 5))

plt.plot(
    temperatures,
    Cv,
    "o-"
)

plt.xlabel("Temperature T")
plt.ylabel(r"$C_V/N$")
plt.title(f"Heat capacity per spin, last {N_LAST} frames")
plt.tight_layout()
plt.savefig("cv_vs_temperature.png", dpi=200)
plt.show()


# ------------------------------------------------------------
# Mean acceptance fraction vs T
# ------------------------------------------------------------

plt.figure(figsize=(7, 5))

plt.plot(
    temperatures,
    mean_acceptance,
    "o-"
)

plt.xlabel("Temperature T")
plt.ylabel("Mean acceptance fraction")
plt.title(f"Mean acceptance fraction, last {N_LAST} frames")
plt.tight_layout()
plt.savefig("acceptance_vs_temperature.png", dpi=200)
plt.show()

