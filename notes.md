### General definitions
If we can define a phase, then the variables that describe the kinematics of a plane wave might be written as a function of the phase $S$:

$$
S = {\bf k}\cdot {\bf x} - \sigma t,
$$

where $\bf{k}$ is the wavenumber vector, which is everywhere perpendicular to lines of constant phase, ${\bf x}$ is the position, and $\sigma$ is the wave frequency.  The wavenumber vector and wave frequency are related to the phase by:

$$
{\bf k} = {\bf \nabla} S,
$$
$$
\sigma = - \frac{\partial S}{\partial t}.
$$

it follows that the wavenumber vector is **irrotational**

$$
{\bf \nabla} \times {\bf k} = 0.
$$

Additionally, differentiating the equations for the wavenumber and frequency in time and space, respectively, we gets:

$$
\frac{\partial {\bf k}}{\partial t} + {\bf \nabla}\omega =0,
$$

which is the conservation of wave crests. If there are no waves being created by a local disturbance, the local rate of change in wavenumber (number of waves per unit length) is balanced by the convergence of the number of crests passing a fixed point per unit time (frequency).


Now, up to this point, we have not yet considered any physics on our problem. The physics is contained in the dispersion relationship: 

$$
\sigma = \Omega(\bf{k}).
$$

Our goal is to explore the kinematics of ocean surface gravity waves, for which:

$$
\sigma = [gk\tanh{(kH)}]^{1/2}
$$

###   Waves in Inhomogeneous Media
More specifically, we would like to solve the equations describing the kinematics of surface waves propagating through **inhomogeneous** media. If assume that the phase is a slowly varying function of $\bf{x}$ and $t$, then the **local** frequency and wavenumber are related to the phase in the same way as above. Now, granted that the medium is slowly varying, we can assert the dispersion relation still holds locally:

$$
\sigma = \Omega(\bf{k};\bf{x}, t).
$$

Defining the group velocity as 

$$
c_{g_i} \overset{\mathrm{def}}{=} \frac{\partial \sigma}{\partial k_i} = \frac{\partial \Omega}{\partial k_i}, 
$$

and combining with the conservations of wave crests, we obtain

$$
\frac{d \sigma}{dt}  = \frac{\partial \Omega}{\partial t},
$$

$$
\frac{d {\bf k}}{dt}= - {\bf \nabla} \Omega,
$$

where

$$
\frac{d }{dt} \overset{\mathrm{def}}{=} \frac{\partial}{\partial t} + ({\bf c_g \cdot \nabla})
$$

the total derivative as the derivative following the wave group.

### Ray Equations in computer-friendly form (draft, no current)

The wave state is defined by the wavenumber ${\bf k} = (k_x, k_y)$. Then, we can compute:

$$
k = (k_x^2 + k_y^2)^{1/2}
$$
$$
\theta = \arctan{\left(\frac{k_y}{k_x}\right)}
$$

$$
c_g = \frac{g}{2} \frac{\tanh{(kH)}+ kH \text{sech}^2{(kH)}}{\left[gk\tanh{(kH)}\right]^{1/2}}
$$

$$
c_{g_x} = c_g \cos{\theta}
$$

$$
c_{g_y} = c_g \sin{\theta}
$$

$$
\frac{d x}{dt}= c_{g_x}
$$


$$
\frac{d y}{dt}= c_{g_y}
$$

$$
\frac{d k_x}{dt}= - \frac{1}{2}gk[\tanh{(kH)}]^{-1/2}\frac{\partial H}{\partial x},
$$

$$
\frac{d k_y}{dt}= - \frac{1}{2}gk[\tanh{(kH)}]^{-1/2}\frac{\partial H}{\partial y},
$$

# Current effects
$$
c_{g_x} = c_g  \cos(\theta) + u_x
$$
$$
c_{g_y} = c_g  \sin(\theta) + u_y
$$
$$
\frac{dk_x}{dt}=\frac{dk_x}{dt}(bathy) - k_x\frac{\partial u_x}{\partial x} - k_y\frac{\partial u_y}{\partial x}
$$
$$
\frac{dk_y}{dt}=\frac{dk_y}{dt}(bathy) - k_x\frac{\partial u_x}{\partial y} - k_y\frac{\partial u_y}{\partial y}
$$