#include "r_tf.h"
#include "r_pwm.h"

using namespace Rcpp;

RCPP_EXPOSED_CLASS(PWMPtr)



RCPP_MODULE(mod_tf)
{
  class_<TFPtr>("TF")
  
  .constructor()
  
  .property("name",         &TFPtr::get_name        , &TFPtr::set_name         )
  .property("binding_size", &TFPtr::get_binding_size, &TFPtr::set_binding_size ) 
  .property("coefs",        &TFPtr::get_coefs       , &TFPtr::set_coefs        )
  .property("lambda",       &TFPtr::get_lambda      , &TFPtr::set_lambda       )

  .property("pwm", &TFPtr::get_pwm)
  
  
  ;
}


