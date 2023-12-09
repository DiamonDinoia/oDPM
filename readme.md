# `Fast & accurate Dose Planning Method` R&D

## Description

This project aims to provide a Monte Carlo simulator for photon and electron real-time radiotherapy treatment.
The real-time part of the project needs to be able to compute the radiation doses in at most `1 second`.
To achieve this time requirement, this software implements gpu acceleration is used.

This software is <span style="color:red">NOT VALIDATED</span> for clinical use and there is no warranty for this
software.

## Dependencies

The code uses C++17 features, mainly std::variant to avoid virtual methods on the GPU where they can cause issues.
This also allows to perform a minor optimization where different particles can be stored sequentially in memory in the
same queue, the virtual method alternative does not allow that as a pointer to the base class would be stored. This will
cause a jump in memory.

This project is based on [dpm-g4cpp](https://github.com/mnovak42/dpm-g4cpp) and the
following [fork](https://github.com/DiamonDinoia/dpm-g4cpp) is required to
generate the material data.

The CMAKE file will download dpm-g4cpp and execute the dpm_GenerateData to generate the material data.
At the moment there is support for:

1. G4_WATER
2. G4_AIR

### Other libraries utilized are:

- OpenMP for parallelism
- [mdspan](https://github.com/kokkos/mdspan.git) to handle multi dimensional data
- [optional](https://github.com/DiamonDinoia/Optional) to use optional values used to detect if the particle is out of
  bounds in an elegant way i.e:

```cpp 
if(geometry[particle.position]) {
 particle inside the geometry
} else {
 out of bounds
}
``` 

### GPU Dependencies:

- [cuda-api-wrappers](https://github.com/eyalroz/cuda-api-wrappers.git) C++ wrapper of CUDA API, this is used to avoid
  ugly marcos such as
  ```cpp
  #define gpuErrorAssert(ans) code = ans; std::cerr << "GPU ERROR: " << cudaGetErrorString(code) << " " << __FILE__
                  << " " << __LINE__ << std::endl; std::abort(-1);

  ```
  since the library will throw proper c++ exceptions in case of error;
- **CUDA >= 9** NVIDIA GPU programming API
- **cuRAND** Shipped with CUDA

### Optional Dependencies (disabled by default):

- [PCG](https://www.pcg-random.org/index.html) fast and reliable random number generator;
- [XoshiroCpp](https://github.com/Reputeless/Xoshiro-cpp) fast (maybe fastest currently) and reliable random number
  generator (reccomented on the CPU);
- [RPMALLOC](https://github.com/mjansson/rpmalloc) high performance memory allocator. It slightly increases performance.
- [CATCH2](https://github.com/catchorg/Catch2) used to unit test the codebase against g4-dpm.

The cmake file will download all the necessary dependencies.
dpm-g4cpp depnds on Geant4 so make sure that Geant4 is installed and source the necessary files before trying to build
this project

## Build

1. As mentioned before, the data generation and simulation phases are well separated. Since the data generation part
   depends on `Geant4`, a relatively recent version (>= Geant4-10.7p01) of the `Geant4` simulation toolkit needs to be
   available on the system. See
   the **[Geant4 Installation Guide](https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/index.html)**
   for further instruction.
2. Options:
    1. -DCMAKE_BUILD_TYPE=Release : release build with compiler optimizations enabled;
    2. -DCMAKE_BUILD_TYPE=Debug : release build with compiler optimizations disable and debug symbols;
    3. -DCMAKE_BUILD_TYPE=RelWithDebInfo : release build with compiler optimizations enabled and debug symbols;
    4. -DGPU=ON enable GPU acceleration, requires CUDA
    5. -DENABLE_LTO=ON enables link time optimizations, minor performance improvement but greatly increases compilation
       time;
    6. -DUSE_PCG=ON enables the use of PCG random number generator on the CPU
    7. -DUSE_XOSHIRO=ON enables the use of PCG random number generator on the CPU (reccomended)
    8. -DRPMALLOC=ON enables the use of RPMALLOC memory allocator. Mainly used to check memory leaks in developement
       more than for perfmance improvements
    9. -DTESTS=ON compiles the Catch2 unit tests
    10. -DPYTHON_WRAPPER=ON (WIP) enables the compilation of a python wrapper. DO NOT ENABLE, broken at the moment

## USAGE

At the moment a parametric main with predefined configuration is not supported but the code can be used as a library.
The following steps are necessary to do a simulation:

1. Define the size of the geometry:
    ```cpp
    ThreeVector<unsigned long> geometrySize{dimx, dimy, dimz};
    ```
   NOTE: the size is in voxel, this corresponds to mm only id the resolution is {1,1,1}
2. Define the size of the voxel (the resolution of the simulation:
    ```cpp
    ThreeVector<real_type> resolution(voxel_size, voxel_size, voxel_size);
    ```
   NOTE: the resolution is in mm
3. Decide the beam for the simulation between:
    ```
    PhotonPencilBeam, ElectronPencilBeam, ToyPhotonGun and ToyElectronGun
    ```
4. Define the beam parameters
    1. Pencil Beams need initial position, initial direction and initial energy:
       ```cpp
       ThreeVector<real_type> initialPosition{dimx * .5 * voxel_size, dimy * voxel_size * .5, theRZ0};
       ThreeVector<real_type> initialDirection{0, 0, 1};
       const auto primaryEnergy = 6.;
       ```
       NOTE: The vectors are in mm, dimx * voxel_size converts the size of the geometry in mm and *.5 calculates the
       center
       of the geometry and the energy is in Mev.
    2. *ToyBeams* requires radiotherapy-style parameters:
        1. RZ0 distance from the geometry
        2. Energy initial energy of the particles
        3. source to axis distance (sad) in mm
        4. field of view (ffl) in mm
        5. offset center of the beam  
           example (to set it to the center of a VoxelCube geometry):
       ```cpp
           ToyElectronGun beam{theRZ0, primaryEnergy, 1000, 100, {real_type(dimx * .5 * voxel_size),
                 real_type(dimy * .5 * voxel_size)}};
       ```
5. Define the Geometry. Currently, there are two geometries supported. Both of them are rectangular 3D geometries
    1. **HalfDistanceVoxelCube**: Rectangular geometry the center is 0,0 and coordinates are in the
       range [-size/2, size/2]
    2. **VoxelCube**: Rectangular geometry the center is [size/2, size/2] and coordinates are in the range [0, size]
       They are completely equivalent, VoxelCube should be slightly faster as the distance to boundary computation is
       simpler.
6. Instantiate the Run to perform the initialization. Run requires two template parameters `<Beam, Geometry>` which
   should have been selected according to the instruction above.
   `Run<T,V>` executes the simulation on the CPU in parallel. `GPURun<T,V>` instead executes it on the GPU.
   The parameters of Run are:
    1. histories: the number of primary histories to simulate
    2. geometry size: ThreeVector with the size of the geometry in voxel
    3. beam an instance of the beam, it should match the template parameter
    4. resolution: resolution of the simulation in mm
    5. seed: seed of the simulation
    6. data directory: directory containing the material data generated by g4-dpm the cmake file defined DATA_DIR and
       stores data into it for commodity;
    7. [optional] Parallelism number of CPU threads or number of GPUs to use (multiple GPUs are supported).
7. Run the simulation and get the results:
    ```cpp
    auto voxelMap = run.simulate();
    ```
   NOTE: VoxelMap is a 3D array containing the radiation doses. For now they are normalized by doing voxelMap/=max(
   voxelMap). This should be eventually converted to <sup>MeV</sup>&frasl;<sub>mm<sup>3</sup></sub>
8. Get the dose depth distribution, integrate x,y plane over the z-axis
    ```cpp
   const auto result = voxelMap.doseDepthDistribution();
    ```
9. Get the dose depth profile, the values seen at the center [x,y] of the geometry over the z-azis. Useful in
   radiotherapy as the experimental
   data contains only this measurement. This is the vector voxelMap(center_x, center_y). This should be eventually
   converted to <sup>MeV</sup>&frasl;<sub>mm<sup>3</sup></sub>
   ```cpp
   const auto result = voxelMap.doseDepthProfile();
   ```

## TESTING

1. compile with `-DTESTS=ON`
2. `cd [BUILD DIRECTORY]/tests & ctest` (might take hours)

## TODO

- [ ] Add a configuration to the fork that generates all the biological NIST materials, Titanium, Water and Air.
- [ ] Add an API where materials can be specified per Voxel. `simulation()` now sets all the voxels to Water.
- [ ] Change the results from pure numbers to <sup>MeV</sup>&frasl;<sub>mm<sup>3</sup></sub>
- [ ] Implement the python wrapper for easy plotting
- [ ] Parse argc and argv to allow selecting different pre-defined
- [ ] Now data is generated for every build but this is not necessary, find a better way to generate data only if there
  is a change to the materials
- [ ] download and compile Geant4 (?) maybe if requested by the user (they have a github repo with a cmake)

ACKNOWLEDGEMENTS:

1. topas for providing
   [ks-test.cpp](https://github.com/tumi8/topas/blob/master/detectionmodules/statmodules/wkp-module/ks-test.cpp)
