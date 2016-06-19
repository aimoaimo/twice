// Crystallographic Inside-Out Projections
// POVRay/MegaPOV 
// aimo.winkelmann@gmail.com
 
               
// MegaPOV required for user_defined camera               
#version unofficial MegaPov 1.21;    
global_settings { max_trace_level 20 } 
#include "math.inc"    
#include "colors.inc"
#include "textures.inc"    
                         
                           
// apply user-defined projection camera (true) or standard orthographic camera (false)                          
#declare CrystalCamera=true;   

#declare AtomScale =0.42; // scaling of atom size relative to touching sphere radius  (close packing)  
#declare Magnification =1.5;       // 1.0 = upper hemisphere in stereographic   


              
                           
#declare Crystal = "BCC"; 
//#declare Crystal = "FCC";   
//#declare Crystal = "HCP"; 
//#declare Crystal = "HEX";  // primitive hexagonal
//#declare Crystal = "SC";   
//#declare Crystal = "TRK"; 

// Euler Angles in Degrees   
#declare  EulerZ1 = 0;     
#declare  EulerX  = 0;   
#declare  EulerZ2 = 0;  



/*
// HCP 001      
#declare  EulerZ1 = 0.0;     
#declare  EulerX  = 0.0;    
#declare  EulerZ2 = 180.0; 
*/
             
             
// FCC 111  
/*    
#declare  EulerZ1 = 45.0;     
#declare  EulerX  = 54.7;    
#declare  EulerZ2 = 0.0;  
*/
                   
#declare Projection ="STEREO";
//#declare Projection ="GNOMONIC";
//#declare Projection ="LAMBERT";   
//#declare Projection ="ORTHO";  


// explore site symmetry: shift origin (in terms of unit cell lattice vectors)  
#declare ashift=0.0;     
#declare bshift=0.0;
#declare cshift=0.0;


// BCC octahedral interstitial                  
//#declare ashift=0.5;     
//#declare bshift=0.5;
//#declare cshift=0.0;    


/*
// BCC tetrahedral interstitial
#declare ashift=0.00;     
#declare bshift=0.50;
#declare cshift=0.25;
*/

// FCC tetrahedral interstitial
//#declare ashift=0.25;     
//#declare bshift=0.25;
//#declare cshift=0.25;


#declare n_max=15;   // maximum number of +/- unit cells for making of crystal


/*
// appearance of atoms
#declare AtomTexture =  texture { pigment{ 
                                            // color rgb<1,1,1>
                                  CoolCopper // MediumOrchid //  Salmon //GreenYellow // Goldenrod //Salmon //SteelBlue //Wheat //White // SlateBlue //OldGold // BrightGold // CoolCopper //RichBlue //MediumOrchid// BrightGold // OrangeRed 
                                  }
                            finish {  phong 0.4 phong_size 100   // phong_size larger for smaller spot
                                     specular 0.3 roughness 0.025      // more rough for larger spot
                                     diffuse 0.35 ambient 0.0
                                     reflection  { 0.45 /*metallic*/ }
                            }
                         } 
*/

// appearance of atoms
#declare AtomTexture =  texture { pigment{ 
                                            // color rgb<1,1,1>
                                            RichBlue // Salmon // CoolCopper //RichBlue //MediumOrchid// BrightGold // OrangeRed 
                                  } 
                                  finish { 
                                            phong 0.4 phong_size 40   // phong_size larger for smaller spot
                                            specular 0.25 roughness 0.5      // more rough for larger spot
                                            diffuse 0.3 ambient 0.15
                                            reflection { 0.25 }
                                            //irid { 0.35 thickness 0.5 turbulence 0.5 } 
                            }
                         } 




        
        
           
// We keep the radius R of the atoms constant for comparison between fcc and bcc
// This means that the lattice constant a has to be scaled to keep the atomic density constant.
// bcc: atoms touch along cube diagonal: 4R=sqrt(3)*a_bcc
// fcc: atoms touch along face diagonal: 4R=sqrt(2)*a_fcc
// which results in: a_fcc=sqrt(3)/sqrt(2)                   
// see: Goldstein "Physical Foundations of Materials Science", Chapter 2            
#declare a_bcc=1.0;
#declare a_fcc=sqrt(3)/sqrt(2);       
              
// TOUCHING-Atom-Radius R (R_fcc = R_bcc due to scaling of a)                                  
#declare R_bcc=0.25*sqrt(3)*a_bcc;                   
#declare R_fcc=0.25*sqrt(2)*a_fcc;   
#declare R_touch=R_fcc; // also =R_bcc

#declare a_hcp=2*R_fcc; // hcp a lattice vector is 2x Touching-Atom-Radius      
#declare c_hcp=sqrt(8/3)*a_hcp;                     
                     
              
              
              
// SETUP of crystal lattice parameters here   
#switch (0)
#case (strcmp(strupr(Crystal),"BCC"))
        // bcc crystal
        #declare a=a_bcc;
        #declare b=a_bcc;
        #declare c=a_bcc;  
        
        #declare uc_alpha=90.0*pi/180;
        #declare uc_beta =90.0*pi/180;
        #declare uc_gamma=90.0*pi/180; 
         
#break  
#case (strcmp(strupr(Crystal),"SC"))
        #declare a=2*R_touch;  
        #declare b=a;
        #declare c=a;
        #declare uc_alpha=90.0*pi/180;
        #declare uc_beta =90.0*pi/180;          
        #declare uc_gamma=90.0*pi/180; 
#break   
#case (strcmp(strupr(Crystal),"TRK"))     // triclinic bravais lattice
        #declare a=2*R_touch;  
        #declare b=a;
        #declare c=a;
        #declare uc_alpha=70.0*pi/180;
        #declare uc_beta =22.0*pi/180;          
        #declare uc_gamma=90.0*pi/180; 
#break 
#case (strcmp(strupr(Crystal),"FCC"))
        // fcc crystal               
        #declare a=a_fcc;
        #declare b=a_fcc;
        #declare c=a_fcc;
        #declare uc_alpha=90.0*pi/180;
        #declare uc_beta =90.0*pi/180;          
        #declare uc_gamma=90.0*pi/180;     
#break 
#case (strcmp(strupr(Crystal),"HCP"))
        // hcp crystal               
        #declare a=a_hcp;
        #declare b=a_hcp;
        #declare c=c_hcp;
        #declare uc_alpha=90.0*pi/180;
        #declare uc_beta =90.0*pi/180;          
        #declare uc_gamma=60.0*pi/180; 
#break
#case (strcmp(strupr(Crystal),"HEX"))
        // hex crystal               
        #declare a=a_hcp;
        #declare b=a_hcp;
        #declare c=c_hcp;
        #declare uc_alpha=90.0*pi/180;
        #declare uc_beta =90.0*pi/180;          
        #declare uc_gamma=60.0*pi/180; 
#break  
#else
      #debug "Not a dog, cat or ostrich.\n"
        // Stuff to do if an unexpected creature was encountered.
#end







              
// O'Keeffe Chapt.4, 4.5, p.109     

// unit cell volume
#declare V_cell=a*b*c*sqrt(1+2*cos(uc_alpha)*cos(uc_beta)*cos(uc_gamma)-cos(uc_alpha)*cos(uc_alpha)-cos(uc_beta)*cos(uc_beta)-cos(uc_gamma)*cos(uc_gamma));
                                         
#declare a_star=b*c*sin(uc_alpha)/V_cell;
#declare b_star=a*c*sin(uc_beta )/V_cell;                                         
#declare c_star=a*b*sin(uc_gamma)/V_cell;
          
#declare cos_alpha_star=(cos(uc_beta)*cos(uc_gamma)-cos(uc_alpha))/(sin(uc_beta)*sin(uc_gamma));          

// Cartesian system p. 112, (4.17)
#declare ax=a;               #declare ay= 0.0;                           #declare az=0.0;
#declare bx=b*cos(uc_gamma); #declare by= b*sin(uc_gamma);               #declare bz=0.0;
#declare cx=c*cos(uc_beta) ; #declare cy=-c*cos_alpha_star*sin(uc_beta); #declare cz=1.0/c_star;






                                 




#declare AtomRadius=AtomScale*R_touch; // size of atoms                                                             
#declare r_outer=12.0*a_fcc; // outer crystal sphere radius in lattice constants


          
           
// define ultility functions for stereographic/gnomonic projection

// u,v -> projection        

#declare s2p = function(s) { Magnification*(-1.0+2.0*s)} 


#declare phi   = function(gx,gy) {     atan2(gy,gx)    }   

// length of vector in projection plane  
#declare rxy    = function(gx,gy) {sqrt(f_sqr(gx)+f_sqr(gy))} 
#declare rxy1   = function(gx,gy) {min(sqrt(f_sqr(gx)+f_sqr(gy)),1.0)} // limit to 1.0
#declare rxy2   = function(gx,gy) {sqrt(f_sqr(gx)+f_sqr(gy))*( max(abs(gx),abs(gy))/rxy(gx,gy) )} // adjust to square as function of Phi
                      
                      
                      
                      
// select type of projection here  
#switch (0)
#case (strcmp(strupr(Projection),"STEREO"))
        // #debug "Assuming STEREOGRAPHIC projection\n"
        #declare theta = function(gx,gy) {2.0* atan(rxy(gx,gy))}   
         
#break 
#case (strcmp(strupr(Projection),"GNOMONIC"))
        // #debug "Assuming GNOMONIC projection\n"
        #declare theta = function(gx,gy) {1.0* atan(rxy(gx,gy))}   
#break
#case (strcmp(strupr(Projection),"LAMBERT"))   
        // Lambert azimuthal equal-area projection of sphere  
        // #debug "Assuming LAMBERT projection\n"
        #declare theta = function(gx,gy){ 2.0*acos(rxy1(gx,gy)) }    
#break  
#else
      #debug "Not a dog, cat or ostrich.\n"
        // Stuff to do if an unexpected creature was encountered.
#end  

  

                         
// x,y,z unit vector components from theta (polar distance)  and phi (azimuth)                         
#declare xpol  = function(gx,gy) {sin(theta(gx,gy))*cos(phi(gx,gy))  }
#declare ypol  = function(gx,gy) {sin(theta(gx,gy))*sin(phi(gx,gy))  }
#declare zpol  = function(gx,gy) {cos(theta(gx,gy))                  }
        

// set up camera type
#if (CrystalCamera)
        camera{
          user_defined
          location{
            function{0}  // location of camera always in origin 
            function{0}  
            function{0}   
          }
          direction{
            function{xpol(s2p(u),s2p(v))}
            function{ypol(s2p(u),s2p(v))}
            function{zpol(s2p(u),s2p(v))}
          }                  
        }
#else
        camera {
                perspective 
                //orthographic
                location <0,0,0>  
                up    <0,1,0>
                right  <1,0,0>
                look_at <0,0,10>
                angle 90
                }  
#end    









#macro PutAtom (afrac,bfrac,cfrac)      
  
  // find atom position           
  #declare Px=(afrac+ashift)*ax+(bfrac+bshift)*bx+(cfrac+cshift)*cx;
  #declare Py=(afrac+ashift)*ay+(bfrac+bshift)*by+(cfrac+cshift)*cy;
  #declare Pz=(afrac+ashift)*az+(bfrac+bshift)*bz+(cfrac+cshift)*cz;
  
  
  // distance of current lattice point from origin                 
  #local r_current=sqrt(f_sqr(Px)+f_sqr(Py)+f_sqr(Pz)); 
      
  #if (r_current>AtomRadius & r_current<r_outer) 
     object{
        sphere {<Px,Py,Pz>,AtomRadius 
                        material{ texture { AtomTexture }} 
        }    
            
        // ZXZ Euler Angle Rotation sequence    
        rotate  EulerZ1  * z     
        rotate  EulerX   * x    
        rotate  EulerZ2  * z  
                      
     } 
  #end          
#end // ------------------ end of macro
       
          
          
  
          
        
#declare ic=-n_max;
#while(ic <= n_max)         
   
  #declare ib = -n_max; // place unit cells +/- n_max times 
  #while(ib <= n_max)

   
   #declare ia = -n_max; 
   #while(ia <= n_max)

        #switch (0)
        #case (strcmp(strupr(Crystal),"BCC"))
                // #debug "Assuming BCC crystal\n"
                // Settings to be used for bcc
                PutAtom(ia,ib,ic)
                PutAtom(ia+0.5,ib+0.5,ic+0.5)  
                 
        #break  
        #case (strcmp(strupr(Crystal),"SC"))
                //#debug "Simple Cubic\n"
                // Settings to be used for simple cubic
                PutAtom(ia,ib,ic)
        #break 
                #case (strcmp(strupr(Crystal),"SC"))
                //#debug "Triclinic\n"
                // Settings to be used for triclinic
                PutAtom(ia,ib,ic)
        #break 
        #case (strcmp(strupr(Crystal),"FCC"))
                //#debug "Face Centered Cubic\n"
                // Settings to be used for FCC
                PutAtom(ia    ,ib    ,ic    )
                PutAtom(ia+0.5,ib+0.5,ic+0.0)
                PutAtom(ia+0.0,ib+0.5,ic+0.5)
                PutAtom(ia+0.5,ib+0.0,ic+0.5)     
        #break 
        #case (strcmp(strupr(Crystal),"HCP"))
                //#debug "HCP cell\n"
                // Settings to be used for HCP
                PutAtom(ia    ,ib    ,ic    )
                PutAtom(ia+1/3,ib+1/3,ic+0.5)
        #break
        #case (strcmp(strupr(Crystal),"HEX"))
                //#debug "primitive HEX cell\n"
                // Settings to be used for HEX
                PutAtom(ia    ,ib    ,ic    )
        #break  
        #else
              #debug "Not a dog, cat or ostrich.\n"
                // Stuff to do if an unexpected creature was encountered.
        #end   
      
      #declare ia=ia+1;
   #end  // ia   
   
   #declare ib=ib+1;
  #end // ib

#declare ic=ic+1;
#end //ic
                
                 
 

// Light source at origin 
light_source {<0,0,0>  White
    fade_distance 4
    fade_power 2
    //shadowless
} 


// add fog for making atoms lighter with distance
/*
fog {
    distance 5  // fade distance of fog
    color rgb<1.0, 1.0, 1.0>
    turbulence 0.5
    turb_depth 3.4
  } 
*/
  
 
/* 
// other lights
// Light source at origin 
light_source {<0,0,0> LightWood // MandarinOrange SlateBlue LightWood OrangeRed Flesh
    fade_distance 1.5
    fade_power 3
    //shadowless   // shadowless will not show phong und specular
} 


// add fog for making atoms lighter with distance
fog {
    distance 30  // fade distance of fog
    color rgb<1.0, 1.0, 0.7>
    turbulence 23.5  // 23.5
    turb_depth 10.4  // 10.4
  }
*/  

