#include <Rcpp.h>
#include "r_defaults.h"

using namespace Rcpp;

Defaults::Defaults()
{
  cex      = 1.25;
  cex_axis = 1.25;
  cex_main = 1.5;
  cex_lab  = 1.5;
  cex_sub  = 1.25;
  lwd      = 2;
  tcl      = 0.5;
  xaxs     = 'i';
  yaxs     = 'i';
}

RCPP_MODULE(mod_par_defaults)
{
  class_<Defaults>("par_defaults")
  
  .constructor()
                                                                      
  .property("cex"     , &Defaults::get_cex     , &Defaults::set_cex     )  
  .property("cex.axis", &Defaults::get_cex_axis, &Defaults::set_cex_axis) 
  .property("cex.main", &Defaults::get_cex_main, &Defaults::set_cex_main)
  .property("cex.lab" , &Defaults::get_cex_lab , &Defaults::set_cex_lab )
  .property("cex.sub" , &Defaults::get_cex_sub , &Defaults::set_cex_sub )
  .property("lwd"     , &Defaults::get_lwd     , &Defaults::set_lwd     )
  .property("tcl"     , &Defaults::get_tcl     , &Defaults::set_tcl     )
  .property("xaxs"    , &Defaults::get_xaxs    , &Defaults::set_xaxs    )
  .property("yaxs"    , &Defaults::get_yaxs    , &Defaults::set_yaxs    )
  ;
}
