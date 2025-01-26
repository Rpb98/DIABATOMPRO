# DIABATOMPRO
<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a id="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
<!-- [![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url] -->
<!-- [![MIT License][license-shield]][license-url] -->
<!-- [![LinkedIn][linkedin-shield]][linkedin-url] -->



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="images/logo.png" alt="Logo" width="750" height="500">
  </a>

  <!-- <h3 align="center">Best-README-Template</h3> -->

  <p align="center">
    A simple to use Julia package to jumpstart your diabatisation project!
    <br />
    <a href="#quick-example-usage"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <!-- <a href="https://github.com/othneildrew/Best-README-Template">View Demo</a> -->
    <!-- · -->
    <a href="https://github.com/othneildrew/Best-README-Template/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/othneildrew/Best-README-Template/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installing-julia">Installation</a></li>
        <li><a href="#quick-example-usage">Quick Example Usage</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <ul>
        <li><a href="#the-input-file">The Input File</a></li>
            <ul>
                <li><a href="#the-bond-grid">The Diatomic Bond Grid</a></li>
                <li><a href="#the-method-block">The Diabatisation Method Block</a></li>
                <li><a href="#the-save-block">Saving Your Diabatisation</a></li>
                <li><a href="#the-potential-block">Potential Energy Curves</a></li>
                <li><a href="#the-nac-block">Non-Adiabatic Couplings (NACs)</a></li>
                <li><a href="#the-dipole-block">Dipoles</a></li>
                <li><a href="#the-spin-orbit-block">Spin-Orbit Couplings</a></li>
                <li><a href="#the-electronic-angular-momentum-block">Electronic Angula Momentum</a></li>
                </ul>
          <li><a href="#fitting-2-state-system-nacs">Fitting NACs: The 2-State Approximation</a></li>
    </ul>
    <li><a href="#contact">Contact Details</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>

                


  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
[![Product Name Screen Shot][CH-diab]](https://www.exomol.com/)

The (stationary) Schrodinger equation for atomistic systems is solved using the adiabatic potential energy curves (PECs) and the associated adiabatic approximation. Despite being very simplistic, this approach is very powerful and used in nearly all practical applications where the equillibrium properties of molecules are usually well represented. In cases when interactions between electronic states become important, the associated  non-adiabatic effects are taken into account via the derivative couplings (DDRs), also known as non-adiabatic couplings (NACs). For diatomic molecules, the corresponding PECs in the adiabatic representation are characterized by avoided crossings. The alternative to the adiabatic approach is the diabatic representation, obtained via a unitary transformation of the adiabatic states by minimizing the DDRs. For diatomics, the diabatic representation has zero DDR and non-diagonal diabatic couplings (DCs) ensue. The two representations are fully equivalent and so should be the rovibronic energies and wavefunctions which result from the solution of the corresponding Schrodinger equations.

DIABATOMPRO is a julia code that computes the adiabatic to diabatic transformation (AtDT) for the diatomic N-coupled electronic state problem. Property based and hybrid asymptotic-property-based approaches are employed to determine the AtDT, with an easy to use interface for initialising such calculations. Furthermore, a regualarisation scheme for the non-adiabatic couplings is implemented, where they are usually not internally consistent with eachother or the associated molecular properties such as potentials and dipoles. As a result, a non-physical diabatic representation of the system ensues. Thus, the aim of the regularisation procedure is to ensure physical diabatic properties upon transformation of the adiabatic system by correction of the NACs. For a more detailed description of the problem, please see the articles: [An ab initio study of the rovibronic spectrum of sulphur monoxide (SO): diabatic vs. adiabatic representation](https://pubs.rsc.org/en/content/articlelanding/2022/cp/d2cp03051a), [Numerical Equivalence of Diabatic and Adiabatic Representations in Diatomic Molecules](https://pubs.acs.org/doi/10.1021/acs.jctc.3c01150), [papers in preperation](https://www.researchgate.net/profile/Ryan_Brady13).


<p align="right">(<a href="#readme-top">back to top</a>)</p>


### Dependencies

This package relies on the following Julia libraries:

- `Dierckx`
- `Interpolations`
- `LinearAlgebra` (Standard Library)
- `ApproxFun`
- `Optim`
- `StaticArrays`
- `QuadGK`
- `Distributed` (Standard Library)
- `Statistics` (Standard Library)
- `SpecialFunctions`
- `Trapz`
- `DataFrames`
- `CSV`

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

Before using this package, you need to have the following:

- **Julia**: Version 1.6 or higher is recommended.
- A working knowledge of the **Julia package manager**.

### Installing Julia

To download and install Julia, follow these steps:

1. Go to the official Julia website: [https://julialang.org/downloads/](https://julialang.org/downloads/)
2. Choose the installer for your operating system (Windows, macOS, or Linux).
3. Follow the instructions on the website to install Julia on your system.
4. Once installed, ensure Julia is added to your system path by running `julia` in your terminal or command prompt.

### Installing DIABATOMPRO

Once Julia is installed, you can install and use the `DIABATOMPRO` package by following these steps:

1. Open a terminal or command prompt.
2. Start Julia by typing `julia` and pressing enter.
3. In the Julia REPL, use the package manager to add the package:

    ```julia
    ] add https://github.com/Rpb98/DIABATOMPRO.git
    ```

4. Once installed, you can start using the package in your Julia scripts by adding:

    ```julia
    using DIABATOMPRO
    ```

### Quick example Usage

After installing the package, here's an example of how to run the main function:

```julia
using DIABATOMPRO

# Example of using the main function
U, Adiabatic_Objects, Diabatic_Objects, Diabatic_Basis, Hamiltonian = Diabatise("path/to/input-file.inp", save_flag = true , special_name = "test01")
```
After running the above function, a diabatic representation of your input adiabatic system should be computed. Please see a comprehensive guide below for more details. The function has the following inpit & output syntax:

**Input:**
* `"path/to/input-file.inp"` : the path to the input file ([see here for details on the input file](#the-input-file))

* `save_flag` : if set to **true** then the diabatisation will be saved according to the parameters in the input file ([see here for details on saving](#the-save-block))

* `special_name` : is a string that will be incoorporated into the saved output files of the diabatisation.

**Output:**
* `U` : the adiabatic to diabatic transformation (AtDT) with shape 

* `Adiabatic_Objects` : a dictionary of adiabatic property matrices, each of shape 

* `Diabatic_Objects` : a dictionary of diabatic property matrices, each of shape 

* `Diabatic_Basis` : a vector of size $N$, each element corresponding to the r-dependent basis coefficients in the adiabatic basis, 

```math
\begin{align*}
\ket{\psi^{\rm(d)}_i}=\sum^N_j\mathcal{C}_{ij}\ket{\psi^{\rm(a)}_j}
\end{align*}
```

* `Hamiltonian` : a dictionary of Hamiltonian elements, each a structs which holds information about the adiabatic objects defined in the input file.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage
### The Input File
The input file for **DIABATOMPRO** follows closesly the input file for the rovibronic code **Duo** ([for more details please see the following documentation/article](https://github.com/ExoMol/Duo)). Blocks are defined for the calculation parameters, save file parameters, and molecular property curves in the *adiabatic representation*.

#### The Bond Grid
The grid of bond lengths can be specified via the following example syntax:

```
grid                
  npoints 2001
  range   0.7 3.7   
  interpolation_type  cubic
end
```
* `npoints` - is the number of equally spaced points defining the grid which properties are interpolated onto.
* `range` - the range (in any unit) of the bond length grid, in this case the diatom bond lengths start at 0.7 and end at 3.7 with 2001 equally spaced points in between.
* `interpolation_type` - currently only cubic has been implemented. Quintic interpolation is currently under development.

#### The Method Block
The calculations meta-variables can be defined in the following block:
```
method   
    atoms 14N 14N             
    states 1 2 3
	diabatisation forward-evolution
	regularisation potential
	grid_resolution 0.0001
	r_boundary_condition 0 0 -1 1 0 0 0 -1 0
end
```
* `atoms` - the reduced mass of the molecule is computed based on the atomic masses of the input atoms. In this example the 14N (natural) isotope is specified, which a mass of 14.003074 amu will be retrieved from `atomic_masses.jl`. Or, one can specify their own mass directly as a float.

* `states` - the state numbers of the states wanted to be included in the diabatisation procedures. These numbers are usually assigned to electronic states on energy enumeration of the potentials, but can be any number the user wants (i.e. can skip intermediate states).

* `diabatisation` - the diabatisation algorithms currently implemented are:
    * `2-state-approx`: When dealing with only two coupled electronic states (approximately isolated states), an analytical adiabatic to diabatic unitary transformation can be generated ([see here for details](https://pubs.acs.org/doi/10.1021/acs.jctc.3c01150)). This means only two states in the `states` card can be specified.

    * `forward-evolution`: An AtDT is computed via the formal exponential line-integral propagator method. The `forward` card means a solution for the AtDT is evolved forwards from a boundary condition at short bond lengths (usually the identity) to larger stretches. If the user does not want to evolve from the identity at short bond lengths, instead they can specify a boundary condition by adding the line `l_boundary_condition` (see card later for syntax).

    * `backward-evolution`: An AtDT is computed via the formal exponential line-integral propagator method. The `backward` card means a solution for the AtDT is evolved backwards from a boundary condition at long bond lengths (usually a signed permutation matrix) to shorter stretches. The user can specify the large-stretch boundary condition by the card `r_boundary_condition` (see card later for syntax). If this is left blank, instead a forwards evolved solution is computed, and the closest signed permutation matrix to the forwards solution at the boundary is taken.

    * `evolution`: Two AtDTs are computed via the formal exponential line-integral propagator method via a `forward`  and `backward` evolution from the physical (or specified) boundary conditions. Because of the inconsistency of the NACs between eachother and the molecular properties, the `forward`  and `backward` evolved solutions do not neccesarily align. To fix this, a regualrising correction to the NACs are computed via connection of the associated generator matrices (Lie algebras $\mathfrak{so}(N)$) such that the diabatic properties remain smooth. As a result, a smooth set of diabatic properties are ensured to have the correct assymptotic limits, a new set of corrected NACs are computed, and a new AtDT which satisfies both boundary conditions is generateed. [A papers detailing the method is currently in preperation](https://www.researchgate.net/profile/Ryan_Brady13).
    
* `regularisation` - the molecular property used in the regularisation procedure provoked by the `diabatisation evolution` card who's diabatic representation will be made smoot. This can be any object of the following: `potential`, `dipole`, `spin-orbit`, `lx`. 

* `grid_resolution` - parameter specifying the coursness of the grid used in solution of the AtDT through evolution. It gives the order of magnitude in grid point seperation. This grid usually needs to be very dense for an accurate solution to the AtDT ($\sim10^{-5}$ Å) and would be inefficient/cumbersome to define a grid in the `grid` block with such a small grid seperation. The AtDT is then splined onto the global grid defined in the `grid` block.

* `r_boundary_condition` / `l_boundary_condition` - the matrix elements of the desired boundary values of the AtDT. **r** means the long bond length (or 'right') boundary and **l** means the short bond length (or 'left') boundary. The matrix elements are listed with space delimiters, counting left/right and top/down, for example a $3\times3$ matrix would have elements listed as $M_{11}$, $M_{12}$, $M_{13}$, $M_{21}$, $M_{22}$, $M_{23}$, $M_{31}$, $M_{32}$, $M_{33}$. e.g. `0 0 -1 1 0 0 0 -1 0` as in the example would yield the following matrix:
```math
\begin{align*}
\left(\begin{matrix}0 & 0 & -1 \\ 1 & 0 & 0 \\ 0 & -1 & 0\end{matrix} \right)
\end{align*}
```
#### The Save Block
The results of the diabatisation can be saved in different formats, and are specified via the following block:
```
save
	as duo
	jrot 0
	vmax 1000
	nroots 250
end
```
* `as` - defines the output format of the diabatisation results and can be set as:

    * `duo`: saves diabatic and adiabatic representations of the diatomic spectroscopic model in the form of a **Duo** input file ([see here for details & documentation](https://github.com/ExoMol/Duo)). In this case, extra parameters can be set:

        * `jrot`: the possible values of the total angular momentum quantum number $J$ for the **Duo** calculation.

        * `vmax`: the contracted vibrational basis size to be used in the **Duo** rovibronic calculations.

        * `nroots`: the number of eigenvalues (rovibronic energy terms) of the diatomic Schrodinger equation of molecular motion to compute by **Duo**.

    * `DataFrame`: saves all of the adiabatic & diabatic properties as `.csv` files in the current directory.

#### The Potential Block
Adiabatic potential energy curves can be defined by the following block:

```
poten 1
name (1)^1Sigma+
symmetry +
lambda 0
mult 1
type grid
units bohr cm-1
values
    1.499999999282	226784.954900000000
    1.504999999179	223868.253900000000
    .               .
    .               .
    .               .        
end
```
* `poten` - this is the potential keyword used to define a potential object. Following this is a unique integer identifier for the state, and typically follows the energy enumeration of the electornic states input to the program. However, one can use any integer here to label their state. Note: this number is directly the row/column index of the state in the Hamiltonian matrix, i.e. `poten 1` would place the defined potential at index $V^{(\rm a)}_{11}$ in the adiabatic potential matrix. This identifier will be used in all subsequent couplings.

* `name` - the user defined name for the potential, it will be used as the name in the saved **Duo** input file if saved.

* `symmetry` - the parity of the electronic state, it defines whether the electronic wavefunction changes phase upon inversion of the electron coordinates through the molecular center of mass. If $+$, then the electron distribution is symmetric about this point, if $-$ then it is antisymmetric. It is used in this code to determine whether it should have an avoided crossing with another nearby state.

* `lambda` - the projected total orbital angular momenta onto the internuclear axis of the diatomic. It is used in this code to determine whether it should have an avoided crossing with another nearby state.

* `mult` - the spin multiplicity, defined as $2S+1$. It is used in this code to determine whether it should have an avoided crossing with another nearby state.

* `type` - is the functional type of the input object. Currently, grid representation is suitable for defining complex adiabatic potentials, and so one should use `type`. In principle one can use functional forms, but this feature and its uses (fitting & modelling adiabatic potentials using diabatic potentials and couplings) are in development. 

* `units` - the units of the bond length then potential. All curves are converted to Å and $\rm cm^{-1}$. See `units` on the possible conversions currently available.

* `values` - begins the grid/function parameter block. After this, space delimited columns represent the bond length vs. potential.

* `end` - marks the end of the `poten` block.

The electronic term symbols `symmetry`, `lambda`, and `mult` are checked against other states to ensure diabatisation is performed on states of _the same symmetry_ **only**.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

#### The NAC Block
The Non-Adiabatic Couplings (NACs) define the diabatisation procedures. and can be defined by the following blocks.

**As a grid**:

```
NAC i j
name <i|d/dr|j>
spin Si Sj
sigma Σi Σj
lambda Λi Λj
type <type>
sub-types <type1> <type2> <type3>
factor 1.0
values 
    .           .
    .           .
    .           .
end
```
* `NAC` - this is the non-adiabatic coupling keyword used to define a NAC object. Following this are two integer ID's for the bra and ket state, an are linked to the states $i$ and $j$ of the `poten` objects. These numbers are the row/column index of the bra and ket state in the Hamiltonian matrix, i.e. `NAC 1 2` would place the defined NAC at index $W^{(\rm 1)}_{12}$ in the NAC (first derivative coupling) matrix.

* `name` - the user defined name for the NAC, it will be used as the name in the saved **Duo** input file if saved.

* `spin` - the spin of the bra and ket electronic states (should be the same).

* `sigma` - the projection of the total spin angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `lambda` - the projection of the total orbital angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `factor` - a scalar factor applied to the NAC coupling curve.

* `sub-types` - the functional type of the NAC sometimes requires subtypes. Please see `mix` and `mix_perturbed` below for examples.

* `type` - the functional type of the NAC input to the program. It can take the following entries:

`grid` - grid representation of the NAC, given by two columns of bond length (Å) versus NAC (1/Å) in the `values` block (see potentials for an example of a grid input).

`lorentzian` - a lorentzian profile is used to model the NAC, and is often used for its cusp-like shape. It is parameterised by the Half Width at Half Maximum (HWHM), $\gamma$, the peak position, $r_0$, and amplitude $N$.  It is programmed as:

  <!-- $$f(r;\gamma,r_0)=\frac{N}{2}\frac{\gamma}{\gamma^2+(r-r_0)^2},$$ -->
<p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728735700.jpg" alt="equation" />
</p>


and enters the input file via the following:
```
NAC 1 2
name <1|d/dr|2>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type lorentzian
factor 1.0
values 
    gamma 0.01466346 
    r0    1.22717420
    N     1.02726695
end
```

>_NOTE: if modelling a 2-state problem, N=1 should be used as this ensures the correct mixing at the crossing point and correct asymptotic limits_

`laplacian` - a laplacian profile is used to model the NAC, and is desireable for its cusp-like shape. Typically it overestimates the NAC at the peak and underestimates the NAC in the wings. It is parameterised by the Half Width at Half Maximum (HWHM), $\gamma$, the peak position, $r_0$, and amplitude $N$.  It is programmed as:

 <!-- $$f(r;\gamma,r_0)=\frac{N\pi}{4\gamma}\exp\left(-\frac{|r-r_0|}{\gamma}\right),$$ -->

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728735778.jpg" alt="equation" />
</p>

and enters the input file via the following:
```
NAC 1 2
name <1|d/dr|2>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type laplacian
factor 1.0
values 
    gamma 0.01466346 
    r0    1.22717420
    N     1.02726695
end
```

>_NOTE: if modelling a 2-state problem, N=1 should be used as this ensures the correct mixing at the crossing point and correct asymptotic limits_

`gaussian` - a gaussian profile is used to model the NAC. It is parameterised by the width, $\gamma$, the peak position, $r_0$, and amplitude $N$. It is programmed as:

<!-- $$f(r;\gamma,r_0)=\frac{N}{2\gamma}\exp\left(-\ln(2)\left(\frac{r-r_0}{\gamma}\right)^2\right),$$ -->


 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728735816.jpg" alt="equation" />
</p>

and enters the input file via the following:
```
NAC 1 2
name <1|d/dr|2>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type gaussian
factor 1.0
values 
    gamma 0.01466346 
    r0    1.22717420
    N     1.02726695
end
```

`lor_lap` - the geometric average of a lorentzian and laplacian mixing angle (cumulative distribution function) is used to model the NAC ([see here for details](https://pubs.rsc.org/en/content/articlelanding/2022/cp/d2cp03051a)). It is parameterised by two widths, $\gamma$ & $\delta$, the peak position, $r_0$, and is normalised to unit area. It is programmed as:

<!-- $$f(r;\gamma,r_0)=\frac{d\beta^{\rm avg}_{ij}}{dr}$$ -->

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751210.jpg" alt="equation" />
</p>

where the mixing angle between states $i$ and $j$, $\beta^{\rm avg}_{ij}$, is

<!-- $$\beta^{\rm avg}_{ij} = \frac{1}{2}\arcsin\left(\sqrt{\sin(2\beta^{\rm lo}_{ij})\sin(2\beta^{\rm la}_{ij})}\right)$$ -->

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751280.jpg" alt="equation" />
</p>

where the lorentzian and laplacian mixing angles are

<!-- $$\beta^{\rm lo}_{ij}=\frac{\pi}{4}+\frac{1}{2}\arctan\left(\frac{r-r_0}{\gamma}\right)$$
$$\beta^{\rm la}_{ij}=\begin{cases} 
\frac{\pi}{4}\exp(\frac{r-r_0}{\delta}) & \text{if } r < r_0, \\
\frac{\pi}{4} & \text{if } r = r_0, \\
\frac{\pi}{2}-\frac{\pi}{4}\exp(-\frac{r-r_0}{\delta}) & \text{if } r > r_0.
\end{cases}$$ -->

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751389.jpg" alt="equation" />
</p>

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751415.jpg" alt="equation" />
</p>
where the maximal overlap of the lorentzian and laplacian function is ensured when

$$\delta = 1.397\times \gamma$$

and is used for this function. Therefore only $\gamma$ is required as input by the user, with the following syntax:

```
NAC 1 2
name <1|d/dr|2>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type lor_lap
factor 1.0
values 
    gamma 0.01466346 
    r0    1.22717420
end
```

`mix` - a linear combination of two functions is used to model the NAC. An additional line in the NAC block must be used to specify which two functions are mixed with the following syntax:
```
sub-types f1 f2
```
In this case the two functions $f_1$ and $f_2$ are mixed in the following linear combination

<!-- $$f = mf_2 + (1-m)f_1$$ -->

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751530.jpg" alt="equation" />
</p>

where $m$ controls the fraction of the final function character contributed by either sub-functions. $f_1$ and $f_2$ are allowed to be any programmed function, and share the same peak position $r_0$, width $\gamma$, and amplitude $N$. The  `mix` function allows for greater flexibilty by including assymetry into the function profile by the allowing the width parameter to vary sigmoidally as

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751599.jpg" alt="equation" />
</p>

<!-- $$\gamma(r;a,r_0)=\frac{2\gamma_0}{1+\exp(a(r-r_0))}$$ -->

meaning at the crossing point the width parameter reduces to a reference width $\gamma_0$. $a>0$ will give the function a long tail left of the peak, $a<0$ will give the function a long tail right of the peak, and $a=0$ makes the function symmetric about $r_0$. 

An example input block would look like:

```
NAC 2 3
name <2|d/dr|3>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type mix
sub-types lorentzian gaussian
factor 1.0
values 
    gamma  0.027265880
    r0     1.115074610
    N     -1.011149990
    a     -11.81517805
    m      0.392749990
end
```


`mix_perturbed` - a linear combination of two functions with sigmoidally varying width parameter (skewness) and an arbitrary set of perturbing functions is used to model the NAC. An additional line in the NAC block must be used to specify which two functions are mixed and what functional type of perturbing functions with the following syntax:
```
sub-types f1 f2 f3
```
In this case the two functions $f_1$ and $f_2$ are mixed in the following linear combination and perturbation by $f_3$

<!-- $$f(r;r_0,\tilde{\gamma},N) = mf_2(r;r_0,\tilde{\gamma},N)+ (1-m)f_1(r;r_0,\tilde{\gamma},N)+\sum^n_if_3(r;r_{0,i},\gamma_{0,i},N_{0,i})$$ -->

 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751637.jpg" alt="equation" />
</p>

where $m$ controls the fraction of the final function character contributed by either sub-functions. $f_1$, $f_2$, and $f_3$ are allowed to be any programmed function, and share the same peak position and width. The  `mix` function allows for greater flexibilty by including assymetry into the function profile by the allowing the width parameter to vary sigmoidally as

<!-- $$\tilde{\gamma}(r;a,r_0)=\frac{2\gamma_0}{1+\exp(a(r-r_0))},$$ -->
 <p align="center">
  <img src="https://www.sciweavers.org/download/Tex2Img_1728751599.jpg" alt="equation" />
</p>

meaning at the crossing point the width parameter reduces to a reference width $\gamma_0$. $a>0$ will give the function a long tail left of the peak, $a<0$ will give the function a long tail right of the peak, and $a=0$ makes the function symmetric about $r_0$. 

An example input block would look like:

```
NAC 1 3
name <1|d/dr|3>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type mix_perturbed
sub-types lorentzian gaussian lorentzian
factor 1.0
values 
    gamma   0.21776621
    r0      1.30726219
    N       0.16291701
    a     -18.71219532
    m      -0.59683866
    gp      0.0300       (first perturbing function)
    r0p     1.0300                   .
    Np      0.0130                   .
    gp      0.0150       (second perturbing function)
    r0p     1.0970                   .
    Np      0.0065                   .
    gp      0.0150       (third perturbing function)
    r0p     1.1400                   .
    Np      0.0035                   .
    .       .                        .
    .       .                        .
    .       .                        .
end
```
#### The Dipole Block
The adiabatic dipoles are input similarly to the other objects in **DIABATOM-PRO**, and can be defined by the following block:
```
dipole i j
name <i|DMC|j>
spin Si Sj
sigma Σi Σj
lambda Λi Λj
lz     lz_i lz_j
type <type>
units angstrom Debye
factor 1.0
values 
    .           .
    .           .
    .           .
end
```
* `dipole` - this is the dipole moment coupling keyword used to define a dipole object. Following this are two integer ID's for the bra and ket state, an are linked to the states $i$ and $j$ of the `poten` objects. These numbers are the row/column index of the bra and ket state in the Hamiltonian matrix, i.e. `dipole 1 2` would place the defined NAC at index $\mu_{12}$ in the dipole matrix.

* `name` - the user defined name for the dipole, it will be used as the name in the saved **Duo** input file if saved.

* `spin` - the spin of the bra and ket electronic states (should be the same).

* `sigma` - the projection of the total spin angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `lambda` - the projection of the total orbital angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `lz` - the elements of the $L_z$ operator, used in **Duo** to convert cartesian representation to $\Lambda$ representation. It is only used when aaving diabatisation results to a duo input file.

* `type` - the functional type of the dipole input to the program. It can take the following entries:

`grid` : grid representation of the dipole, given by two columns of bond length (Å) versus dipole (Debye) in the `values` block (see [potentials](#the-potential-block) for an example of a grid input).

<!-- `polynom_decay_24` : -->

#### The Spin-Orbit Block
The adiabatic spin-orbit couplings are input similarly to the other objects in **DIABATOM-PRO**, and can be defined by the following block:
```
spin-orbit i j
name <i|DMC|j>
spin Si Sj
sigma Σi Σj
lambda Λi Λj
lz     lz_i lz_j
type <type>
units angstrom cm-1
factor 1.0
values 
    .           .
    .           .
    .           .
end
```
* `spin-orbit` - this is the spin-orbit coupling keyword used to define a spin-orbit coupling object. Following this are two integer ID's for the bra and ket state, an are linked to the states $i$ and $j$ of the `poten` objects. These numbers are the row/column index of the bra and ket state in the Hamiltonian matrix, i.e. `spin-orbit 1 2` would place the defined spin-orbit at index $SO_{12}$ in the spin-orbit matrix.

* `name` - the user defined name for the spin-orbit coupling, it will be used as the name in the saved **Duo** input file if saved.

* `spin` - the spin of the bra and ket electronic states (should be the same).

* `sigma` - the projection of the total spin angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `lambda` - the projection of the total orbital angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `lz` - the elements of the $L_z$ operator, used in **Duo** to convert cartesian representation (like in **MOLPRO** outputs) to $\Lambda$ representation. It is only used when aaving diabatisation results to a duo input file.

* `type` - the functional type of the spin-orbit input to the program. It can take the following entries:

`grid` : grid representation of the spin-orbit coupling, given by two columns of bond length (Å) versus spin-orbit ($\rm cm^{-1}$) in the `values` block (see [potentials](#the-potential-block) for an example of a grid input).

<!-- `polynom_decay_24` : -->

#### The Electronic Angular Momentum Block
The adiabatic spin-orbit couplings are input similarly to the other objects in **DIABATOM-PRO**, and can be defined by the following block:
```
lx i j
name <i|DMC|j>
spin Si Sj
lambda Λi Λj
lz     lz_i lz_j
type <type>
units angstrom cm-1
factor 1.0
values 
    .           .
    .           .
    .           .
end
```
* `lx` - this is the electronic angular momentum (EAM) coupling $L_x$ keyword used to define an EAM object. Following this are two integer ID's for the bra and ket state, an are linked to the states $i$ and $j$ of the `poten` objects. These numbers are the row/column index of the bra and ket state in the Hamiltonian matrix, i.e. `lx 1 2` would place the defined EAM at index $L_{12}$ in the EAM matrix.

* `name` - the user defined name for the EAM coupling, it will be used as the name in the saved **Duo** input file if saved.

* `spin` - the spin of the bra and ket electronic states (should be the same).

* `sigma` - the projection of the total spin angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `lambda` - the projection of the total orbital angular momentum on the internuclear axis of the bra and ket electronic states (should be the same).

* `lz` - the elements of the $L_z$ operator, used in **Duo** to convert cartesian representation (like in **MOLPRO**) to $\Lambda$ representation. It is only used when saving diabatisation results to a duo input file.

* `type` - the functional type of the EAM input to the program. It can take the following entries:

`grid` : grid representation of the EAM, given by two columns of bond length (Å) versus EAM (units of $\hbar$) in the `values` block (see [potentials](#the-potential-block) for an example of a grid input).

<!-- `polynom_decay_24` : -->


### Fitting 2-State System NACs
In the case when one deals with only two-states that exhibit an avoided crossing, one can fit the NACs to ensure a smooth diabatic representation, known as _property based diabatisation_. First the `diabatisation` [card](#the-method-block) needs to be set to **2-state-approx** to activate the relevent algorithm. Next, if one defines a NAC between states $1$ and $2$ (please see guidlines on defining a NAC [here](#the-nac-block)), for example
```
NAC 1 2
name <1|d/dr|2>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type gaussian
factor 1.0
values 
    gamma 0.01466346 
    r0    1.22717420
    N     1.02726695
end
```
Then a diabatising transformation will be computed by integration of the above NAC to obtain the mixing angle. However, to fit this NAC such that the resulting diabatic potentials are smooth, one can use the `fit` keyword as so
```
NAC 1 2
name <1|d/dr|2>
spin 0.0 0.0
sigma 0.0 0.0
lambda 0 0
type gaussian
factor 1.0
values 
    gamma 0.01466346 fit
    r0    1.22717420 fit
    N     1.00000000
end
```
This will fit only the width $\gamma$ and the peak position $r_0$ of the lorentzian. `fit` can be placed on any parameter, which will be fitted to ensure smooth diabatic potentials.



<!-- ROADMAP -->
<!-- ## Roadmap

- [x] Add Changelog
- [x] Add back to top links
- [ ] Add Additional Templates w/ Examples
- [ ] Add "components" document to easily copy & paste sections of the readme
- [ ] Multi-language Support
    - [ ] Chinese
    - [ ] Spanish

See the [open issues](https://github.com/othneildrew/Best-README-Template/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>
 -->


<!-- CONTRIBUTING
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Top contributors:

<a href="https://github.com/othneildrew/Best-README-Template/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=othneildrew/Best-README-Template" alt="contrib.rocks image" />
</a>

<p align="right">(<a href="#readme-top">back to top</a>)</p>
 -->


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Ryan Brady - ryan.brady.17@ucl.ac.uk


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

[![Product Name Screen Shot][group]](https://www.exomol.com/)

My name is Ryan, I am a third yeard PhD student at the ExoMol group (UCL). ExoMol is a database of molecular line lists that can be used for spectral characterisation and simulation - and as input to atmospheric models of exoplanets - of brown dwarfs and cool stars, and other models including those for combustion and sunspots. The ExoMol format and database is described in the recent article, [The 2024 release of the ExoMol database: Molecular line lists for exoplanet and other hot atmospheres](https://www.sciencedirect.com/science/article/pii/S0022407324001900). All the code used in the ExoMol group can be found at the [github repository here](https://github.com/ExoMol).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[CH-diab]: images/CH_diabatisation_adi-dia.png
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
[group]: images/exomol_logo.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 
