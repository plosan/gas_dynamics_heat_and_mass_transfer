# gas_dynamics_heat_and_mass_transfer

This is a repository to store the Gas Dynamics and Heat and Mass Transfer course project.

## The course
Gas Dynamics is a 6th quarter course in the Aerospace Engineering Technologies degree. In it we study the various forms of Heat Transfer, namely,
conduction, convection and radiation. The equations that arise in this field are usually PDEs with difficult and mixed boundary conditions, therefore
we cannot expect to obtain an analytical solution. In order to solve Heat Transfer problem, we follow a numerical approach. We discretize the domain
where the heat transfer problem takes place is discretized and then apply the necessary physical principles (conservation of mass, conservation of
energy, conservation of momentum) to derive the discretization equations. Finally we write an algorithm that solves the problem.

Since the course is focused on the numerical treatment of Heat Transfer, we are encouraged to do an individual project on it.

## The preparation-project
The aim of the preparation-project is to prove one can solve a numerical Heat Transfer problem. Hence an easy case of Heat Transfer that takes short time to
be solved is recommended. In this case a long tube (like a pipe) was considered. The inner and outer faces of the tube are in contact with fluids at
different temperatures which vary over time according to a given function, so it is a transient heat transfer problem. The purpose of the project is
to find the wall temperature as a function of time.

The program to solve the problem was written in MATLAB due to development speed reasons.

## The project
Once the preparation-project is finished and checked by the professors, they propose several problems of heat transfer among which one must choose.

### Problem statement
A 2-dimensional heat conduction problem is considered. The domain is the cross-section of a long bar, so that heat conduction does not occur along
the bar, but in the cross-section. The bar is made of 4 different materials, each one with its own thermophysical properties. Each side of the
cross-section has a different boundary condition.
- On the lower side a constant temperature is imposed (non homogeneous Dirichlet condition).
- On the upper side there is a heat inflow which sets the value of the derivative at the boundary (non homogeneous Neumann condition).
- On the left side a fluid transfers heat according to Newton's law (Robin condition).
- On the right side the temperature increases as a linear function of time (non homogeneous Dirichlet condition).
There are no internal heat sources.

### Numerical aspects
Some numerical aspects of the resolution:
- The domain was discretized using a node-centered mesh with 198 x 200 nodes.
- The time integration was solved with a Crank-Nicolson scheme, as it is a 2nd order unconditionally convergent scheme. The implicit scheme, which
is a 1st order unconditionally convergent scheme is also available.
- At each time instant a linear system has to be solved. The unknowns of this system are the temperatures of each node. A direct method to solve the
problem (such as the LU factorization) was not suitable, since the matrix was sparse. Moreover, the matrix was diagonally dominant by rows. This makes
Gauss-Seidel algorithm ideal to solve the linear system, as the convergence is guaranteed.

### Verification
A verification process is needed in order to assure that the code produces results with physical sense. Two verification processes were carried out.
- The first verification consisted on checking the conservation of energy. This was made separately with an implicit and a Crank-Nicolson scheme, each with
three different time steps. The results were good, with errors of the order of 10^{-10} watts.
- The second verification consisted on producing a plot that showed the isotherm curves and comparing it with the plot given by the professors. This
one was harder. The program had to search the mesh for those nodes with a given temperature. The final plot resembled to that given by the professors,
but it was not perfect. A better algorithm to find the isotherms is the following. The mesh can be thought of as a graph. Hence, when the first node
with the given temperature is found, the algorithm performs a DFS starting at the initial node. Then it visits the neighbor nodes. Since the temperature
must be a continuous function of (x,y), the algorithm goes to the neighbor node with the closest temperature. This enhanced algorithm to find the isotherms
was not implemented due to time reasons.

### Discretization influence analysis
Once the verification is finished some numerical aspects were evaluated. More precisely, the mesh density, time step and time integration scheme were
modified in order to evaluate their influence in the solution. No significant differences between the solutions were found.

### New heat transfer problem
Finally a new heat transfer case was considered. The new domain was the squared cross-section of a long bar, made of a single material. Each side of
the bar was in contact with the same fluid. The initial temperature of the cross-section was 200ºC, while the fluid temperature was 20ºC. The heat
transfer coefficient of the fluid was 100 W/m^2 K.
A radial distribution of temperature varying over time was obtained as expected. As time goes to infinity, the temperature of the cross-section slowly
tends to that of the fluid.

### Code
The code for the project was written in MATLAB. Other languages such as C++ or Fortran are better suited for these kinds of projects due to their
execution speed. Nonetheless, writing code in those languages is slower than in MATLAB, since many features are not implemented. The program is not
scalable, in the sense that it works for the proposed case. In spite of this, the generalization to other cases of 2D heat transfer is pretty
straightforward.

## Some mathematical aspects
As mentioned before, the equations appearing in this course are mainly Partial Differential Equations. In order to solve a heat transfer problem,
initial conditions and boundary conditions are required. Hence, an important part of the subject is solving initial value problems (IVPs) involving the
diffusion equation.

When dealing with IVPs one cares about the existence and uniqueness of a solution. As stated above, the equation is the diffusion equation
<img src="https://render.githubusercontent.com/render/math?math=u_t + D \Delta u = 0"> in the domain
<img src="https://render.githubusercontent.com/render/math?math=\Omega = (0,W) \times (0, H) \subset \mathbb{R}^2"> with mixed boundary conditions at
<img src="https://render.githubusercontent.com/render/math?math=\partial \Omega">.

Consider the general homogeneous diffusion problem at <img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}^n \times [0, +\infty)">
with initial condition <img src="https://render.githubusercontent.com/render/math?math=g in C(\mathbb{R}^n) \cap L^\infty(\mahtbb{R}^n)">. A well-known
result is that the solution to the problem is a <img src="https://render.githubusercontent.com/render/math?math=g in C^\infty(\mathbb{R}^n \times (0,\infty))">
function (Theorem 1, section 2.3, Partial Differential Equations - Lawrence C. Evans).

By doing reflections and periodic extensions, it might be possible to extend the 2D heat transfer problem to a problem in
<img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}^2 \times [0, +\infty)">. This extended problem might not be homogeneous,
that is, it might have internal heat sources given by a function <img src="https://render.githubusercontent.com/render/math?math=f \colon \mathbb{R}^2 \times [0,+\infty) \rightarrow \mathbb{R}">.
If <img src="https://render.githubusercontent.com/render/math?math=f"> is good enough, then the extended problem has existence and uniqueness of solution.
The restriction of this solution to <img src="https://render.githubusercontent.com/render/math?math=\Omega"> would be the solution to the initial problem.
