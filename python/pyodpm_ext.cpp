#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include "DPMG4Run.h"
#include "beams.h"
#include "run.h"

namespace nb = nanobind;
using namespace nb::literals;

using OPMCG4Run = opmc::G4Run;

NB_MODULE(pyodpm_ext, m) {
    m.attr("G4DATASET_DIR") = G4DATASET_DIR;
    m.attr("MeV")           = CLHEP::MeV;
    m.doc()                 = "A simple example python extension";
    nb::class_<opmc::ThreeVector<double>>(m, "ThreeVector_d")
        .def(nb::init<double, double, double>())
        .def_rw("x", &opmc::ThreeVector<double>::x)
        .def_rw("y", &opmc::ThreeVector<double>::y)
        .def_rw("z", &opmc::ThreeVector<double>::z);
    nb::class_<opmc::ThreeVector<std::int64_t>>(m, "ThreeVector_i")
        .def(nb::init<int64_t, int64_t, int64_t>())
        .def_rw("x", &opmc::ThreeVector<int64_t>::x)
        .def_rw("y", &opmc::ThreeVector<int64_t>::y)
        .def_rw("z", &opmc::ThreeVector<int64_t>::z);
    nb::class_<OPMCG4Run>(m, "G4Run")
        .def(nb::init<std::array<G4long, 2>, const std::string &>())
        .def(nb::init<std::array<G4long, 2>, const std::string &, const std::string &>())
        .def("Run", &OPMCG4Run::Run)
        .def("initializeDicom", &OPMCG4Run::initializeDicom)
        .def("initializeWater", &OPMCG4Run::initializeWater)
        .def("convertGeometry", &OPMCG4Run::convertGeometry)
        .def("center", &OPMCG4Run::center);
    nb::class_<opmc::VoxelCube>(m, "VoxelCube")
        .def("doseDepthDistribution", &opmc::VoxelCube::doseDepthDistribution)
        .def("doseDistribution", &opmc::VoxelCube::doseDistribution)
        .def("resolution", &opmc::VoxelCube::resolution);
    nb::class_<opmc::ElectronPencilBeam>(m, "ElectronPencilBeam")
        .def(nb::init<opmc::ThreeVector<double>, opmc::ThreeVector<double>, double>());
    nb::class_<opmc::PhotonPencilBeam>(m, "PhotonPencilBeam")
        .def(nb::init<opmc::ThreeVector<double>, opmc::ThreeVector<double>, double>());
    nb::class_<opmc::Run<opmc::VoxelCube>>(m, "Run").def(nb::init<std::uint64_t>());
    m.def("simulate", [](opmc::Run<opmc::VoxelCube> &run, const OPMCG4Run &g4Run, opmc::VoxelCube &voxelMap,
                         opmc::ElectronPencilBeam &beam, const std::int32_t histories) {
        const auto tables      = g4Run.getTables();
        auto referenceMaterial = tables->referenceMaterial();
        auto electronData      = tables->getElectronData();
        auto photonData        = tables->getPhotonData();
        run.simulate(voxelMap, beam, referenceMaterial, electronData, photonData, histories);
    });
    m.def("simulate", [](opmc::Run<opmc::VoxelCube> &run, const OPMCG4Run &g4Run, opmc::VoxelCube &voxelMap,
                         opmc::PhotonPencilBeam &beam, const std::int32_t histories) {
        const auto tables      = g4Run.getTables();
        auto referenceMaterial = tables->referenceMaterial();
        auto electronData      = tables->getElectronData();
        auto photonData        = tables->getPhotonData();
        run.simulate(voxelMap, beam, referenceMaterial, electronData, photonData, histories);
    });
    m.def("validate", [](opmc::VoxelCube &a, opmc::VoxelCube &b) {
        const auto gamma = opmc::validateDistributions(a.doseDistribution(), b.doseDistribution());
        return std::make_tuple(gamma.maxRelativeError, gamma.averageRelativeError, gamma.gammaPassingRate,
                               gamma.minGamma);
    });
    m.def("toNdArray", [](opmc::VoxelCube &a) {
        size_t shape[3] = {a.dimension().x, a.dimension().y, a.dimension().z};
        return nb::ndarray<nb::numpy, const double, nb::shape<3, nb::any>>(a.origin_dose(), 3, shape);
    });
}
