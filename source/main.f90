program main
   use blueprint_module

   type(Blueprint) :: details

   details = read_blueprint("D:\Google Drive\Research\Forward ANN\data\test_blueprint.json")

   print *, details%n_input
   print *, details%output_names
   print *, details%input_traits%std
   print *, details%hidden_traits(1)%activation
   print *, details%output_traits(3)%transformation
end program main
