program main
   use test_activation_module
   use test_bivariate_normal_module
   use test_model_module
   use test_shash_module

   call test_activation()
   call test_bivariate_normal()
   call test_model()
   call test_shash()
end program main
