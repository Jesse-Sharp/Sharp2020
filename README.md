# Sharp2019
Matlab implementation of the forward backward sweep method for multiple optimal controls. This code is used to generate results in Figure S6 of the Supplementary material by Sharp et al.

The pay-off considered in this code corresponds to the case where the chemotherapy control is bang-bang and the stem cell transplant control is bounded continuous. 

The supplied parameters in this code specifically generate Figure S6(a). Results in subsequent rows of Figure S6 can be obtained by adjusting parameters a1, a2 and a3 (Line 35 - Line 37). Results in the final column can be obtained by setting Vupper = 0.3 (Line 43) and adjusting a1, a2 and a3 accordingly. Results in the centre column can be obtained by setting Vupper = Inf (Line 43). 

This particular code example is provided as it demonstrates implementation of both bang-bang control and bounded continuous control. To solve for alternative chemotherapy and stem cell transplant control formulations, Line 103 and Line 107 can be modified, respectively. For example, the results in Figure S4 of the Supplementary material by Sharp et al. can be recovered by setting: 

Line 103 to Uupdate = max(min(Uupper,0.5\*(y(1,:).\*kappa.\*Lambda(1,:)+y(2,:).\*Lambda(2,:))/a1),Ulower); 

and Line 107 to Vupdate = max(Vupper\*((sign((a2+Lambda(1,:)))-1)/(-2)),Vlower\*((sign((a2+Lambda(1,:)))+1)/(2))). 
