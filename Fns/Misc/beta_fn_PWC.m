function beta = beta_fn_PWC(parameter_vector)

%---------------------------------------

liquid_density = parameter_vector(1);
gravity = parameter_vector(5);

%-------------------------------------

beta = Flex_PWC(parameter_vector) / (liquid_density * gravity) ;