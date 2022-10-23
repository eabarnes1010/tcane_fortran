program main
   use test_bivariate_normal_module
   use test_model_module
   use test_shash_module

   write(*,*)
   call test_shash()
   write(*,*)
   call test_bivariate_normal()
   write(*,*)
   call test_model()
end program main
