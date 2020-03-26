#include <Rcpp.h>

#ifndef R_DEFAULTS_H
#define R_DEFAULTS_H

class Defaults
{
private:
  double cex;
  double cex_axis;
  double cex_main;
  double cex_lab;
  double cex_sub;
  double lwd;
  double tcl;
  char   xaxs;
  char   yaxs;
  
public:
  Defaults();

  // getters
  double get_cex()      { return cex;      }   
  double get_cex_axis() { return cex_axis; }  
  double get_cex_main() { return cex_main; }
  double get_cex_lab()  { return cex_lab;  }
  double get_cex_sub()  { return cex_sub;  }
  double get_lwd()      { return lwd;      }
  double get_tcl()      { return tcl;      }  
  char   get_xaxs()     { return xaxs;     }
  char   get_yaxs()     { return yaxs;     }

  // setters                           
  void set_cex(double cex)           { this->cex      = cex;      }      
  void set_cex_axis(double cex_axis) { this->cex_axis = cex_axis; } 
  void set_cex_main(double cex_main) { this->cex_main = cex_main; } 
  void set_cex_lab(double cex_lab)   { this->cex_lab  = cex_lab;  }
  void set_cex_sub(double cex_sub)   { this->cex_sub  = cex_sub;  }
  void set_lwd(double lwd)           { this->lwd      = lwd;      }
  void set_tcl(double tcl)           { this->tcl      = tcl;      }
  void set_xaxs(char xaxs)           { this->xaxs     = xaxs;     }
  void set_yaxs(char yaxs)           { this->yaxs     = yaxs;     }
};

  
  
  
  
  
  
  
  
  
#endif
