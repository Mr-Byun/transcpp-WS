#include <Rcpp.h>
#include "r_mode.h"

using namespace Rcpp;

ModePtr::ModePtr() :
  mode(boost::shared_ptr<Mode>(new Mode())) 
{}


ModePtr::ModePtr(string fname) :
  mode(boost::shared_ptr<Mode>(new Mode())) 
{
  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& input_node = root_node.get_child("Output");

  mode->read(input_node);
}

ModePtr::ModePtr(string fname, string section) :
  mode(boost::shared_ptr<Mode>(new Mode())) 
{
  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& input_node = root_node.get_child(section);

  mode->read(input_node);
}

  
  
RCPP_MODULE(mod_mode)
{
  class_<ModePtr>("Mode")
  
  .constructor()
  .constructor<string>()
  .constructor<string, string>()
  
  .property("occupancy_method", &ModePtr::get_occupancy_method , &ModePtr::set_occupancy_method )    
  .property("score_function"  , &ModePtr::get_score_function   , &ModePtr::set_score_function   )
  .property("scale_data_type" , &ModePtr::get_scale_data_type  , &ModePtr::set_scale_data_type  )
  .property("scale_to"        , &ModePtr::get_scale_to         , &ModePtr::set_scale_to         )
  .property("min_data"        , &ModePtr::get_min_data         , &ModePtr::set_min_data         )
  .property("p_thresh"        , &ModePtr::get_p_thresh         , &ModePtr::set_p_thresh         )
  .property("gc"              , &ModePtr::get_gc               , &ModePtr::set_gc               )
  .property("per_gene"        , &ModePtr::get_per_gene         , &ModePtr::set_per_gene         )
  .property("per_nuc"         , &ModePtr::get_per_nuc          , &ModePtr::set_per_nuc          )
  .property("profiling"       , &ModePtr::get_profiling        , &ModePtr::set_profiling        )
  .property("competition"     , &ModePtr::get_competition      , &ModePtr::set_competition      )
  .property("self_competition", &ModePtr::get_self_competition , &ModePtr::set_self_competition )
  .property("scale_data"      , &ModePtr::get_scale_data       , &ModePtr::set_scale_data       )
  .property("window"          , &ModePtr::get_window           , &ModePtr::set_window           )
  .property("shift"           , &ModePtr::get_shift            , &ModePtr::set_shift            )
  .property("num_threads"     , &ModePtr::get_num_threads      , &ModePtr::set_num_threads      )
  .property("schedule"        , &ModePtr::get_schedule         , &ModePtr::set_schedule         )
  .property("verbose"         , &ModePtr::get_verbose          , &ModePtr::set_verbose          )
  .property("precision"       , &ModePtr::get_precision        , &ModePtr::set_precision        )
  .property("seed"            , &ModePtr::get_seed             , &ModePtr::set_seed             )
  .property("n"               , &ModePtr::get_n                , &ModePtr::set_n                )
  .property("occupancy_method", &ModePtr::get_occupancy_method , &ModePtr::set_occupancy_method )
      
      
  ;
}
  
  
