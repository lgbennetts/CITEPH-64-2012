function w0 = weight_0_PWC(parameter_vector,  root)

h = parameter_vector(7);

w0 =  sech(root*h);