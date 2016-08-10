% AOloopControl
% Olivier Guyon
% Aug 9, 2016

# Overview

## Scope

AO loop control package



## Usage

Scripts to run the software are located within the source code directory:

	./src/AOloopControl/scripts/

The scripts can be linked to your working directory by executing the following command:

	ln -s $PWD/syncscripts /myworkdirectory/syncscripts

Then, execute in your work directory:

	./syncscripts
	
This will install all required scripts in workdirectory and install any packages required.

The main script is

	./aolconf



## Supporting scripts, aolconfscripts directory

Scripts in the `aolconfscripts` directory are part of the high-level ASCII control GUI

------------------------------ -----------------------------------------------------------
Script                         Description
------------------------------ -----------------------------------------------------------
**aolconf_DMfuncs**            DM functions 

**aolconf_DMturb**             DM turbulence functions

**aolconf_funcs**              Misc functions

**aolconf_logfuncs**           data and command logging

**aolconf_menuconfigureloop**  configure loop menu

**aolconf_menucontrolloop**    control loop menu

**aolconf_menucontrolmatrix**  control matrix menu

**aolconf_menu_mkFModes**      Make modes

**aolconf_menurecord**         

**aolconf_menutestmode**       Test mode menu

**aolconf_menutop**            Top level menu

**aolconf_menuview**           Data view menu

**aolconf_readconf**           Configuration read functions

**aolconf_template**           Template (not used)
------------------------------ -----------------------------------------------------------


## Supporting scripts, auxscripts directory

Scripts in the `auxscripts` directory are called by aolconf to perform various tasks.

------------------------------ -----------------------------------------------------------
Script                         Description
------------------------------ -----------------------------------------------------------
**acquRespM**                  Acquire response matrix
------------------------------ -----------------------------------------------------------


## Hardware simulation scripts

SCripts in the `aohardsim` directory are called to simulate hardware for testing / simulations


------------------------------ -----------------------------------------------------------
Script                         Description
------------------------------ -----------------------------------------------------------
**simDMstart**                 Start simulation DM shared mem 
------------------------------ -----------------------------------------------------------




# Simulation environment - PyWFS-based example

## Overview

The AOsim simulation architecture relies on individual processes that simulate subsystems. Each process is launched by a bash script. ASCII configuration files are read by each process. Data I/O can be done with low latency using shared memory and semaphores: a process operation (for example, the wavefront sensor process computing WFS signals) is typically triggered by a semaphore contained in the shared memory wavefront stream. A low-speed file system based alternative to shared memory and semaphores is also provided.


## Processes and scripts

### AO loop control

The ``aolconf`` script is used to configure and launch the AO control loop. It can be configured with input/output from real hardware or a simulation of real hardware.

### Hardware simulation architecture

Close-loop simulation of a pyramid-based WFS requires the following scripts to be launched to simulate the hardware, in the following order :

* ``simDMstart``: This script creates DM channels (uses dm index 5 for simulation). Shared memory arrays ``dm05disp00`` to ``dm05disp07`` are created, along with the total displacement ``dm05disp``. Also creates the ``wf1opd`` shared memory stream which is needed by runDM and will be updated by runWF. ``wf1opd`` is the master clock for the whole simulation, as it triggers DM shape computation and WFS image computation.
* ``AOsim_runDM``: Simulates physical deformable mirror (DM)
* ``AOsim_runWF``: Creates atmospheric wavefronts
* ``AOsim_runPyrWFS``: Simulates pyramid WFS

Some key script variables need to coordinated between scripts. The following WF array size should match :

* ``WFsize`` in script ``simDMstart``
* ``ARRAYSIZE`` in ``AOsim_runWF``
* ``ARRAYSIZE`` in ``AOsim_runDM``


The main hardware loop is between ``AOsim_runWF`` and ``AOsim_runPyWFS``: computation of a wavefront by ``AOsim_runWF`` is *triggered* by completion of a WFS instantaneous image computation by ``AOsim_runPyrWFS``. The configuration files are configured for this link.

The interdependencies between processes in the low-latency shared memory mode are :

~~~
    WFSinst [5] -----(AOsim_runWF)------> wfopd0--->wfopd1=wfopd0-2*dm05dispmap
    dm05dispmap [5]-/

    wfopd1   [5] ----(AOsim_runPyWFS)----> WFinst
    wfopd1   [6] ----(AOsim_runDM)-------> dm05dispmap
~~~

### DM temporal response

The DM temporal response is assumed to be such that the distance between the current position $p$ and desired displacement $c$ values is multiplided by coefficient $a<1$ at each time step $dt$. The corresponding step response is :

$c - p((k+1) dt) = (c - p(k dt)) a$

$c - p(k dt) = (c-p0) a^k$

$p(k dt) = 1-a^k$

The corresponding time constant is

$a^{\frac{t0}{dt}} = 0.5$

$\frac{t0}{dt} ln(a) = ln(0.5)$

$ln(a) = ln(0.5) dt/t0$

$a = 0.5^{\frac{dt}{t0}}$
