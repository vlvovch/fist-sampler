## Calculating acceptance parameters

This repository samples 0-5\% Au-Au collisions from RHIC-BES using MUSIC hypersurfaces and calculates acceptance parameters $\alpha$ for

- net-protons ($0.2 < p_T < 2.0$ GeV/$c$, $|y|$ < 0.5) affected by baryon conservation
- net-$\Lambda$'s ($0.0 < p_T < 2.0$ GeV/$c$, $|y|$ < 0.5) affected by baryon or strangeness conservation
- net-kaons ($0.2 < p_T < 1.6$ GeV/$c$, $|y|$ < 0.5)affected by strangeness conservation
- net-charges ($0.2 < p_T < 2.0$ GeV/$c$, $|\eta|$ < 0.5) affected by charge conservation

The parameter $\alpha$ is defined as ratio of measured particles relative to total conserved charge. For instance for net-protons affected by baryon conservation
$$
\alpha_p^B = \frac{\langle p_{\rm acc}\rangle + \langle \bar{p}_{\rm acc}\rangle}{\langle  B_{4\pi}\rangle + \langle \bar{B}_{4\pi}\rangle}.
$$

### Usage

Run pre-defined scripts. For instance, a GCE calculation within ideal HRG, strong + electromagnetic decays, with 100K events can be run from command line with
```bash
# Compile
make BES-alphas-KpLaQ

# Run
bash alphaBESscanGCE.sh 100000
```
which uses the input files from [input folder](../input)

Modify shell scripts and/or input/source files as appropriate.




