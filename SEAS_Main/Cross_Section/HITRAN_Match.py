"""
list of all hitran molecules and its respective molecule number
more detail see https://hitran.org/lbl/


Molecules not included:
34 O: data only availabe from 68 - 158 wn
37 HOBr: data only availabe from 0.155 - 315 wn
48 C2N2: data only availabe from 200 - 307 wn

Not available via hp.fetch(), other methods possible but not investigated. 
Not urgent since these molecules have low priority.
30 SF6: 
35 ClONO2:
42 CF4: 

"""


HITRAN_Match = {"H2O":1,
                "CO2":2,
                "O3":3,
                "N2O":4,
                "CO":5,
                "CH4":6,
                "O2":7,
                "NO":8,
                "SO2":9,
                "NO2":10,
                "NH3":11,
                "HNO3":12, # very time consuming gas
                "OH":13,
                "HF":14,
                "HCl":15,
                "HBr":16,
                "HI":17,
                "ClO":18,
                "OCS":19,
                "H2CO":20,
                "HOCl":21,
                "N2":22,
                "HCN":23,
                "CH3Cl":24,
                "H2O2":25,
                "C2H2":26,
                "C2H6":27,
                "PH3":28,
                "COF2":29,
                "SF6":30, # not available in hapi
                "H2S":31,
                "HCOOH":32,
                "HO2":33,
                "O":34,
                "ClONO2":35, # not available in hapi
                "NOp":36, # this is NO+
                "HOBr":37, # also doesn't work
                "C2H4":38,
                "CH3OH":39,
                "CH3Br":40,
                "CH3CN":41,
                "CF4":42,# not available in hapi
                "C4H2":43,
                "HC3N":44,
                "H2":45,
                "CS":46,
                "SO3":47,
                "C2N2":48,
                "COCl2":49
                }









