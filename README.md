The current DC optimal power flow (OPF) model finds the optimal operating schedules of generating units considering transmission line constraints.
The objective function is defined as the total operating costs of generating units.
The model type is LP and is applied to a 5 bus PJM test system as described in the following reference:
F. Li and R. Bo, "Small test systems for power system economic studies," IEEE PES General Meeting, Minneapolis, MN, 2010, pp. 1-4.
doi: 10.1109/PES.2010.5589973
The proposed model is able to handle large scale set of data for practical power system transmission networks.

Inputs:
Generator costs + technical characteristics (min/max , connection point)
Network characteristics (reactances, demands, line flow limits)
Outputs:
local marginal prices (LMP)
Generators outputs
Line flows
Angle values
--------------------------------------------------------------------------------------------
Contributed by Dr. Alireza Soroudi, University College Dublin, Dublin, Ireland.
email: alireza.soroudi@gmail.com
