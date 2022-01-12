# Contains some useful compartment numbers


#define the sections in the Acker&Antic2008 model
AckerData = {'soma': range(0,1),'axon': range(1,8), 'basal': range(8,488), 'apical': range(488,1181), 'thetalow':-69}
#Proximal and Distal compartments
#distal_Acker_basal =  [30,52,65,94,114,153,206,226,257,300,323,354,378,401,412,421,461] #30,  
distal_Acker_basal = [206, 487, 153, 461, 226, 257, 300, 378, 23, 30, 323, 354, 421, 52, 65, 94, 401, 412, 114]
#proximal_Acker_basal =  [15,53,66,95,115,158,215,227,258,301,324,470,368,402,413,422,42] #15,  
proximal_Acker_basal =  [154, 462, 115, 422, 214, 227, 263, 368, 21, 24, 302, 324, 413, 41, 53, 73, 386, 402, 111]
distal_Acker_apical = [537, 1149, 1169, 1180, 1130, 1112, 1054, 1092, 1015, 1031]
proximal_Acker_apical = [495, 1131, 1156, 1170, 1113, 1093, 1034, 1055, 996, 1016]
distal_Acker_tuft = [780, 753, 822, 972, 989, 670, 672, 710, 712, 860, 871, 901, 925]
proximal_Acker_tuft = [754, 713, 787, 954, 973, 655, 671, 707, 711, 855, 861, 886, 902]

#define the sections in the Branco&Hausser2011 model
BrancoData = {'soma': range(0,1),'axon': range(1,50), 'basal': range(50,248), 'apical': range(248,488), 'thetalow':-72}
#Proximal and Distal compartments
distal_compartments_Branco = [64,  74,  82,  90,  98, 107, 116, 129, 140, 155, 171, 178,194, 203, 207, 227, 237] #[107,64,194,129,81,171,140,90,247,98,227,153,116,237,203,74,85,155,220,185]
distal_compartments_Branco_eff = [74, 82, 98, 107, 129, 140, 155, 171, 203, 227]
distal_compartments_Branco_nonmda = [64, 90, 116, 178, 194, 207, 237]
proximal_compartments_Branco = [101,55,186,120,66,159,131,238,86,141,92,213,229,108,204]

distal_Branco_basal = [171, 178, 129, 64, 98, 140, 185, 194, 203, 207, 220, 227, 237, 247, 107, 116, 153, 155, 74, 82, 85, 90]    #distal_compartments_Branco_eff
proximal_Branco_basal = [156, 172, 120, 55, 92, 131, 184, 186, 196, 204, 220, 221, 229, 238, 100, 108, 151, 154, 69, 75, 84, 86]    #proximal_compartments_Branco
distal_Branco_apical =  [481, 487, 269, 279, 463, 474, 318, 296, 302, 327, 334]
proximal_Branco_apical =  [476, 482, 261, 270, 449, 464, 303, 293, 297, 321, 328]
distal_Branco_tuft =  [425, 436, 447, 380, 384]
proximal_Branco_tuft =  [385, 431, 437, 363, 381]
#

