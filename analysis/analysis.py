import matplotlib.pyplot as plt
import numpy as np

import pyodpm
import pymedphys


histories = int(1_000_000)

def convertArray(geometry):
    array = pyodpm.toNdArray(geometry)
    return np.ravel(array).reshape(array.shape[0], array.shape[1], array.shape[2], order='F')

def getLateralProfile(array):
    x, y = int(array.shape[0] / 2), int(array.shape[1] / 2)
    return list(array[x, y, :])


def getFrontalDistribution(array):
    return np.sum(np.sum(array, axis=1, initial=0), axis=0, initial=0)


def getSideDistribution(array):
    return np.sum(np.sum(array, axis=2, initial=0), axis=0, initial=0)

def getDepthDistribution(array):
    return np.sum(np.sum(array, axis=0, initial=0), axis=0, initial=0)


run = pyodpm.G4Run([1, 2], pyodpm.G4DATASET_DIR, 'FTFP_BERT_Fast')
# run.initializeWater(333, 333, 333, 1)
run.initializeDicom("/home/mbarbone/repos/photonMonteCarlo/lung_ct/1-", 142)
center = run.center()
print(f'center={center.x}, {center.y}, {center.z}')
resolution = run.resolution()
print(f'resolution={resolution.x}, {resolution.y}, {resolution.z}')

x, y, z = center.x, center.y, center.z

x, y = center.x % resolution.x, center.y % resolution.y
x = resolution.x*.5 - x
y = resolution.y*.5 - y
print(f'x={x}')
print(f'y={y}')
# exit(0)
reference = run.Run("gamma", pyodpm.ThreeVector_d(x, y, -center.z),
                    pyodpm.ThreeVector_d(0., 0., 1.), 6 * pyodpm.MeV, histories)

geometry = run.convertGeometry()

dpm_run = pyodpm.Run(42)
beam = pyodpm.PhotonPencilBeam(pyodpm.ThreeVector_d(center.x+x, center.y+y, .5),
                                 pyodpm.ThreeVector_d(0., 0., 1.), 6 * pyodpm.MeV)
pyodpm.simulate(dpm_run, run, geometry, beam, histories)
# maxError, avgError, gammaPassingRate, minGamma = pyodpm.validate(reference, geometry)

# print(maxError, avgError, gammaPassingRate, minGamma)

reference_array = convertArray(reference)
geometry_array = convertArray(geometry)

plt.plot(getDepthDistribution(reference_array), label="Geant4")
plt.plot(getDepthDistribution(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_depth.png")
plt.draw()
plt.show()
plt.clf()
exit(0)

plt.plot(getLateralProfile(reference_array), label="Geant4")
plt.plot(getLateralProfile(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_lateral.png")
plt.draw()
# plt.show()
plt.clf()

plt.plot(getFrontalDistribution(reference_array), label="Geant4")
plt.plot(getFrontalDistribution(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_frontal.png")
plt.draw()
# plt.show()
plt.clf()

plt.plot(getSideDistribution(reference_array), label="Geant4")
plt.plot(getSideDistribution(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_side.png")
plt.draw()
# plt.show()
plt.clf()

# plt.imshow(geometry_array[:, :, 0], label="ODPM", cmap='Blues')
# plt.draw()
# plt.show()
# plt.clf()
# plt.imshow(reference_array[:, :, 0], label="Geant4", cmap='Reds')
# plt.draw()
# plt.show()
# plt.legend()
# plt.savefig("e_scatter.png")
#
# plt.clf()


dimension = pyodpm.getDimension(reference)
resolution = pyodpm.getResolution(reference)

xmin = 0
ymin = 0
zmin = 0

xmax = dimension.x * resolution.x
ymax = dimension.y * resolution.y
zmax = dimension.z * resolution.z

gridx = resolution.x
gridy = resolution.y
gridz = resolution.z

x = np.arange(xmin, xmax, gridx)
y = np.arange(ymin, ymax, gridy)
z = np.arange(zmin, zmax, gridz)

coords = (x, y, z)

gamma_options = {
    'dose_percent_threshold': 1,
    'distance_mm_threshold': gridz,
    'lower_percent_dose_cutoff': 20,
    'interp_fraction': 10,  # Should be 10 or more for more accurate results
    'max_gamma': 2,
    'random_subset': None,
    'local_gamma': True,
    'ram_available': 2**34,  # 16 GB
    'quiet': True
}

gamma = pymedphys.gamma(
    coords, reference_array,
    coords, geometry_array,
    **gamma_options)
valid_gamma = gamma[~np.isnan(gamma)]

# ---------------------------------------------------------------------
# Passing rate determined here
# ---------------------------------------------------------------------

passing = 100 * np.sum(valid_gamma <= 1) / len(valid_gamma)

# ---------------------------------------------------------------------

plt.figure()
plt.title(f'Gamma Histogram | Passing rate = {passing:.2f}%')
plt.xlabel('Gamma')
plt.ylabel('Number of pixels')

plt.hist(valid_gamma, 20)
plt.show()
