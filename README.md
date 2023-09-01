# Model Predictive Control Scheduling
This project implements a Model Predictive Control (MPC) based scheduling system for optimizing the makespan of jobs on machines. The objective is to efficiently schedule incoming jobs on machines to minimize the overall makespan, which represents the total time required to complete all the jobs.

## How It Works
The scheduling system utilizes a Model Predictive Control approach to make real-time decisions on job allocation to machines. At each arrival of a new job, the system dynamically re-schedules the existing jobs on machines to achieve the best possible makespan

## Test
To test the dynamic and offline scheduling for comparison, you need to:
- Create the **G_init**, **G_j** variables _(G_init0, G_j0 for dynamic scheduling, G_init, G_j for offline scheduling)_
- Create the **BigOmega** values before the optimization loop
- If **BigOmega** has more than one value, choose the element to plot in the **graph_gantt** (_sol\_noNoise(i) and solOpt(i)_ in the last lines of the code in the respective _main\_graph_)

 ### License
 This project is licensed under the terms of the GNU General Public License v3.0
