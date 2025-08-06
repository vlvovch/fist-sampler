# FIST sampler

## Description

FIST sampler is an implementation of the Cooper-Frye particlization sampling procedure for heavy-ion collisions based on a thermal-statistical package [**Thermal-FIST**](https://github.com/vlvovch/Thermal-FIST)

The sampler uses a particlization hypersurface produced in hydro simulations (such as using [**MUSIC**](https://github.com/MUSIC-fluid/MUSIC)) and samples hadrons and resonances from it. Optionally, the chain of decays can then be performed.

A distinctive feature of the sampler is the ability to selectively incorporate the exact (the canonical ensemble) global conservation of conserved charges, such as baryon number, electric charge, strangeness, and charm. It also incorporates the effect of baryon repulsive core (excluded volume), in accordance with the EV-HRG model formulated in [https://arxiv.org/abs/1708.02852](https://arxiv.org/abs/1708.02852)

The sampling procedure first determines the numbers of each hadrons species to be sampled in a given event, and then the momentum of each hadrons is sampled independently from other hadrons. This is different from the more conventional Cooper-Frye sampling and allows one to sample the events with the canonical treatment of conserved charges considerably faster. The downside is that the routine requires a considerable amount of memory, in practice about ~5-15 GB for central Au-Au collisions.

The description of the sampling procedure can be found in Appendix B of the paper [https://arxiv.org/abs/2107.00163](https://arxiv.org/abs/2107.00163)

Currently the sampling does not incorporate viscous corrections. 

## Prerequisites

- A compiler with C++11 support. 

- CMake v3.1+.

## Building

For a Linux system:
```bash
# Clone the repository fetching the submodule(s)
git clone --recurse-submodules https://github.com/vlvovch/fist-sampler.git

# Create the build directory and build the project
cd fist-sampler
mkdir build
cd build
cmake ../
make FISTSampler
```

## Usage
```bash
./FISTSampler <input-file>
```

Here `input-file` contains the parameters for the sampling. 

## Examples

See [input/input.AuAu.7.7.C0-5](input/input.AuAu.7.7.C0-5) for a sample input corresponding to 0-5% central Au-Au collisions at 7.7 GeV, which includes the description of all the relevant parameters.
See also [input/input.AuAu.7.7.C0-5.EVHRG](input/input.AuAu.7.7.C0-5.EVHRG) for a sample input which includes the excluded volume effect for (anti)baryons.

The sampler requires an input file with the particlization hypersurface. One of the supported formats is the binary output from MUSIC.
Such hypersurfaces are available [here](https://drive.google.com/drive/folders/1DMml4IXXcilEZaaTpGF2HM_2ICmeydpz?usp=sharing), based on a paper [https://arxiv.org/abs/2003.05852](https://arxiv.org/abs/2003.05852).
To download the MUSIC hypersurfaces corresponding to 0-5% Au-Au collisions at various RHIC-BES energies run e.g.
```bash
bash ../input/get_MUSIC_input_RHICBES.sh 
```
from the build directory. This will take ~1.5GB of space, this space can be reduced by omitting collision energies that are not needed. Note that the script requires the [gdown](https://github.com/wkentaro/gdown) package to be installed. If this does not work, just download the files from Google Drive manually.

Then, to sample Au-Au collisions at 7.7 GeV run
```bash
./FISTSampler ../input/input.AuAu.7.7.C0-5
```
This should generate 1000 collisions with the canonical treatment of baryon number in ascii output file `AuAu.7.7.C0-5.events.dat`

The [input](input) folder contains also other systems.
For instance, [input/input.HADES.AuAu.2.4.C0-5](input/input.HADES.AuAu.2.4.C0-5) and [input/input.HADES.AuAu.2.4.C0-5.EVHRG](input/input.HADES.AuAu.2.4.C0-5.EVHRG) correspond to Siemens-Rasmussen-Hubble model based event generator of 0-5% central 2.4 GeV Au-Au collisions corresponding to the HADES experiment, with and without the excluded volume effect, respectively.
On other hand, [input/input.ALICE.PbPb.2760.C0-5](input/input.ALICE.PbPb.2760.C0-5) and [input/input.ALICE.PbPb.2760.C0-5.EVHRG](input/input.ALICE.PbPb.2760.C0-5.EVHRG) correspond to blast-wave model based event generator of 0-5% 2.76 TeV Pb-Pb collisions.

## Attribution
Publications using this sampler should include a reference to the **Thermal-FIST** package

- V. Vovchenko, H. Stoecker, *Thermal-FIST: A package for heavy-ion collisions and hadronic equation of state*, [Comput. Phys. Commun. **244**, 295 (2019)](https://doi.org/10.1016/j.cpc.2019.06.024) [[arXiv:1901.05249 [nucl-th]](https://arxiv.org/abs/1901.05249)]

as well as the paper
- V. Vovchenko, V. Koch, C. Shen, *Proton number cumulants and correlation functions in Au-Au collisions at 7.7−200 GeV from hydrodynamics*, [Phys. Rev. C **105**, 014904 (2022)](https://doi.org/10.1103/PhysRevC.105.014904) [[arXiv:2107.00163 [nucl-th]](https://arxiv.org/abs/2107.00163)]

describing the sampling procedure. 

In addition, if third-party input is used (such as the particle list or hypersurface of particlization), please cite the relevant literature.

*Copyright (C) 2021-2022 Volodymyr Vovchenko*
