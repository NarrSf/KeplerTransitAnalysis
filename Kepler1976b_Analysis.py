# Kepler-1976 b
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# 1. Configurations & Constants (This datas are changeable)
DATA_DIR = "/Users/narges_sf/Downloads/kepler1976b_fits"  #File path

PERIOD_DAYS         = 4.959319244       # P [days]
T0_BKJD             = 172.2585292       # Transit epoch (BKJD)
RSTAR_RSUN          = 1.082             # Stellar radius [R_sun]
TRANSIT_DURATION_HR = 2.22739           # [hours]
TRANSIT_DEPTH_PPM   = 9802.0            # catalog depth [ppm]

# 2. Read and Detrend Data 
all_phase = []
all_flux_det = []

for fname in os.listdir(DATA_DIR):
    if not fname.endswith("_llc.fits"):
        continue

    fpath = os.path.join(DATA_DIR, fname)
    with fits.open(fpath) as hdul:
        data = hdul[1].data
        time = data["TIME"]
        flux = data["PDCSAP_FLUX"]
        
    mask = np.isfinite(time) & np.isfinite(flux)
    t = time[mask]
    f = flux[mask]

    if len(t) < 10:
        continue

    # Detrending 
    coefs = np.polyfit(t, f, 3)
    trend = np.polyval(coefs, t)
    f_det = f / trend
    f_det /= np.median(f_det)

    # Phase-folding 
    phase = ((t - T0_BKJD) / PERIOD_DAYS) % 1.0
    phase = phase - 0.5 

    all_phase.append(phase)
    all_flux_det.append(f_det)

phase_all = np.concatenate(all_phase)
flux_all  = np.concatenate(all_flux_det)

# 3. Find Transit Center and Wrap Phase 
nbins = 300
bins = np.linspace(-0.5, 0.5, nbins + 1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
binned_flux = np.array([np.median(flux_all[(phase_all >= bins[i]) & (phase_all < bins[i+1])]) for i in range(nbins)])

idx_min = np.nanargmin(binned_flux)
phi_min = bin_centers[idx_min]

phase_centered = ((phase_all - phi_min + 0.5) % 1.0) - 0.5

# 4. Zoomed Binning for Final Plot 
mask_zoom = (np.abs(phase_centered) < 0.15)
ph_z = phase_centered[mask_zoom]
fl_z = flux_all[mask_zoom]

nbins_z = 100
bins_z = np.linspace(-0.15, 0.15, nbins_z + 1)
bin_centers_z = 0.5 * (bins_z[:-1] + bins_z[1:])
smooth_zoom = np.array([np.median(fl_z[(ph_z >= bins_z[i]) & (ph_z < bins_z[i+1])]) for i in range(nbins_z)])

# 5. Box Model Calculations 
depth_frac = TRANSIT_DEPTH_PPM * 1e-6
duration_phase = (TRANSIT_DURATION_HR / 24.0) / PERIOD_DAYS
half_width = duration_phase / 2.0

phi_model = np.linspace(-0.15, 0.15, 500)
flux_model = np.ones_like(phi_model)
flux_model[np.abs(phi_model) < half_width] = 1.0 - depth_frac

# 6. Visualization 
plt.figure(figsize=(12, 5))

# The main plot (Full Phase)
plt.subplot(1, 2, 1)
plt.scatter(phase_centered, flux_all, s=1, alpha=0.3, color='gray')
plt.plot(bin_centers, binned_flux, color='red', lw=1, label='Binned')
plt.title("Full Phase-Folded Curve (Centered)")
plt.xlabel("Phase")
plt.ylabel("Normalized Flux")

# Zoomed plot (Transit Zoom)
plt.subplot(1, 2, 2)
plt.scatter(ph_z, fl_z, s=5, alpha=0.1, color='orange', label='Data')
plt.plot(bin_centers_z, smooth_zoom, color='red', lw=2, label='Binned Data')
plt.plot(phi_model, flux_model, color='green', lw=2, label='Box Model')
plt.axhline(1.0, ls='--', color='blue', alpha=0.5)
plt.xlim(-0.1, 0.1)
plt.ylim(1.0 - 2*depth_frac, 1.0 + 0.5*depth_frac)
plt.title("Exoplanet Transit (Zoomed)")
plt.xlabel("Phase (centered on 0)")
plt.legend()

plt.tight_layout()
plt.show()

# 7. Final Calculations 
RSUN_IN_REARTH = 109.1
Rstar_Rearth   = RSTAR_RSUN * RSUN_IN_REARTH
Rp_Rearth      = Rstar_Rearth * np.sqrt(depth_frac)

print("-" * 30)
print(f"Stellar Radius: {Rstar_Rearth:.2f} R_earth")
print(f"Calculated Planet Radius: {Rp_Rearth:.2f} R_earth")
print("-" * 30)