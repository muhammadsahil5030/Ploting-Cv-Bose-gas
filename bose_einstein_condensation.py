#______________________________________________importing Modules__________________________________________________
from matplotlib import pyplot as plt
from scipy.integrate import quadrature
from scipy.integrate import quad
from scipy.special import gamma
import numpy as np


T=np.linspace(100,1000,10)
z=[0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

#___________________________________________________******_________________________________________

first_list=[]
second_list=[]
third_list=[]
first1=[0.17995771124404975, 0.36693164091639274, 0.5618622075336955, 0.7659469773712185, 0.9807611476355289, 1.2084682845761785, 1.4522278759444098, 1.71711184388994, 2.0127844892140505, 2.37060366314839]
second2=[0.08147834614093563, 0.16957966841945915, 0.26570891309808364, 0.3718170482488426, 0.49074584857335585, 0.626911765931887, 0.7878508137506574, 0.9884788583012429, 1.2679770552539844, 2.05175480096181]
third3=[0.3383600222311749, 0.7347502340414845, 1.2088139949004502, 1.791119918405047, 2.532521790973207, 3.525694845502204, 4.963522613268125, 7.343650043617557, 12.635329916057854, -4.587839088863846]
squared_list=[0.0066387208898621205, 0.028757263941253712, 0.07060122649976497, 0.13824791736848213, 0.2408314878919831, 0.3930183622638371, 0.620708904727573, 0.9770904533085287, 1.6077658126505658, 4.209697763269837]


for element in second2:
    squared_list.append(element**2)
    division_list=[0.01962028742664641, 0.03913882923599049, 0.058405368234985745, 0.07718518226941995, 0.09509552444934165, 0.11147259745556769, 0.1250541103748252, 0.13305242590606944, 0.12724367494411884, -0.9175774654974541]

division = [i / j for i, j in zip(squared_list, third3)]
div_array=np.array(division_list)
first_array=np.array(first1)


#___________________________________________________performing the Integration________________________________________
for i in z:
    def first_function(x,z):
        gz=((x**(3.0/2))*i*np.exp(-x))/(1-i*np.exp(-x))
        return gz
    result, error = quad( first_function, 0, np.inf, args= (z,))
    Gamma= gamma(5.0/2)
    first_Result= Gamma*result
    first_list.append(first_Result)
   
    def second_function(x,z):
        gz=(x**(1.0/2)*i*np.exp(-x))/(1-i*np.exp(-x))
        return gz
    result, error = quad( second_function, 0, np.inf, args= (z,))
    Gamma= gamma(3.0/2)
    second_Result= Gamma*result
    second_list.append(second_Result)

    def third_function(x,z):
        gz3=(x**(-1.0/2)*i*np.exp(-x))/(1-i*np.exp(-x))
        return gz3
    result, error = quad(third_function, 0, np.inf, args= (z,))
    Gamma= gamma(1.0/2)
    third_Result= Gamma*result
    third_list.append(third_Result)
    
#__________________________________________________print values___________________________________________________

print ("z= ", z , 'Done. \n')
print ("First List= ", first_list, 'Done. \n')
print ("Second List= ", second_list, 'Done. \n')
print ("Third List= ", third_list, 'Done. \n')
print ("Multiply two list= ", squared_list, 'Done. \n')
print ("Division", division, 'Done. \n')
#print ("My array =",division_array, 'Done. \n')
print ("My first array", first_array, 'Done. \n')

#____________________________________________________calculating and Plotting Specific Heat________________________

def specific_heatcapacity(T, first_array, div_array):
    Cv= (3/2)*(T**(3/2))*((5/2)*(first_array)-(3/2)*(div_array))
    return Cv
Cv=np.array(specific_heatcapacity(T, first_array, div_array))
print ("Cv=",Cv,'Done. \n')
plt.plot(T, Cv)
plt.show()
#_____________________________________________________End___________________________________________________________
