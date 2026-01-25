*TITLE eps-Constraint Method for Multiobjective Optimization (EPSCM,SEQ=319)
$ontext
The eps-Constraint Method

$offtext

$inlinecom [ ]
$eolcom //
$STitle Example model definitions

Sets
   i       constraints / 1*1 /
   j       decision variables /1*18 /
   k       objective functions / 1*3 /

Parameter dir(k) direction of the objective functions 1 for max and -1 for min
   / 1  1
     2  1
     3  1
   /

Table c(J,K) matrix of objective function coefficients C
$include "C:\Users\Kostas\Documents\JupyterNotebooks\Rentizelas\CISAMS\pyDecision_Example_WorkedExample_Plastics_PromII\AUGMECON2_4_2_gams\z3_18ctt.txt" ;

Table a(J,I) matrix of technological coefficients A
$include "C:\Users\Kostas\Documents\JupyterNotebooks\Rentizelas\CISAMS\pyDecision_Example_WorkedExample_Plastics_PromII\AUGMECON2_4_2_gams\z3_18att.txt" ;

Parameter b(I)  RHS of the constraints b
/
$include "C:\Users\Kostas\Documents\JupyterNotebooks\Rentizelas\CISAMS\pyDecision_Example_WorkedExample_Plastics_PromII\AUGMECON2_4_2_gams\z3_18btt.txt" ;
/
;
Variables
   Z(K)      objective function variables
Binary Variables
   X(J)    decision variables
Equations
   objfun(K)    objective functions
   con(I)    constraints
;

objfun(K).. sum(J,c(J,K)*X(J)) =e= Z(K);
*con(I).. sum(J, a(J,I)*X(J)) =l= b(I);
con(I).. sum(J, a(J,I)*X(J)) =e= b(I);

display c
display a
display b

Model example / all /;

*---------------------------------------------------------------------------
$STitle eps-constraint method

Set k1(k) the first element of k, km1(k) all but the first elements of k;
k1(k)$(ord(k)=1) = yes; km1(k)=yes; km1(k1) = no;
Set kk(k)     active objective function in constraint allobj

Parameter
   rhs(k)     right hand side of the constrained obj functions in eps-constraint
   maxobj(k)  maximum value from the payoff table
   minobj(k)  minimum value from the payoff table
   numk(k) ordinal value of k starting with 1

Scalar
iter   total number of iterations
infeas total number of infeasibilities
elapsed_time elapsed time for payoff and e-sonstraint
start start time
finish finish time

Variables
   a_objval   auxiliary variable for the objective function
   obj        auxiliary variable during the construction of the payoff table
Positive Variables
   sl(k)      slack or surplus variables for the eps-constraints
Equations
   con_obj(k) constrained objective functions
   augm_obj   augmented objective function to avoid weakly efficient solutions
   allobj     all the objective functions in one expression;

con_obj(km1)..   z(km1) - dir(km1)*sl(km1) =e= rhs(km1);

* We optimize the first objective function and put the others as constraints
* the second term is for avoiding weakly efficient points

* objfun=max z1 + 0.001*(s1/r1+0.1 s2/r2+ 0.01*s3/r3+...)
augm_obj..
  sum(k1,dir(k1)*z(k1))+1.0e-3*sum(km1,power(10,-(numk(km1)-1))*sl(km1)/(maxobj(km1)-minobj(km1))) =e= a_objval;

allobj..  sum(kk, dir(kk)*z(kk)) =e= obj;

Model mod_payoff    / objfun, con, allobj / ;
Model mod_epsmethod / objfun, con, con_obj, augm_obj / ;

Parameter
  payoff(k,k)  payoff tables entries;
Alias(k,kp);

option optcr=0.0;
option limrow=0, limcol=0, solprint=off ;
$offlisting;
$offsymxref;
$offsymlist;
$offuelxref;
$offuellist;
file cplexopt /cplex.opt/;
put cplexopt;
put 'threads 1'/;
put 'parallelmode 1'/;
putclose cplexopt;
mod_epsmethod.optfile=1;
option optca=0.;
mod_payoff.optfile=1;
mod_epsmethod.optfile=1;

* Generate payoff table applying lexicographic optimization
loop(kp,
  kk(kp)=yes;
  repeat
    solve mod_payoff using mip maximizing obj;
    payoff(kp,kk) = z.l(kk);
    z.fx(kk) = z.l(kk); // freeze the value of the last objective optimized
    kk(k++1) = kk(k);   // cycle through the objective functions
  until kk(kp); kk(kp) = no;
* release the fixed values of the objective functions for the new iteration
  z.up(k) = inf; z.lo(k) =-inf;
);
if (mod_payoff.modelstat<>1 and mod_payoff.modelstat<>8, abort 'no optimal solution for mod_payoff');

File fx  / C:\Users\Kostas\Documents\JupyterNotebooks\Rentizelas\CISAMS\pyDecision_Example_WorkedExample_Plastics_PromII\AUGMECON2_4_2_gams\z3_3augm2_2025_f_desvar_Plastics_livedemo.out /;

PUT fx ' PAYOFF TABLE'/   ;
loop (kp,
        loop(k, put payoff(kp,k):12:2);
        put /;
     );
put fx /;

*display payoff;
*minobj(k)=smin(kp,payoff(kp,k));
**$ontext
*minobj(k)=0;
minobj('2')=0;
minobj('3')=0;
maxobj(k)=smax(kp,payoff(kp,k));

*$ontext
*$set fname h.%scrext.dat%

*gridpoints=max integer of km1 = 4149
$if not set gridpoints $set gridpoints 100
Set g grid points /g0*g%gridpoints%/
    grid(k,g) grid
Parameter
    gridrhs(k,g) rhs of eps-constraint at grid point
    maxg(k) maximum point in grid for objective
    posg(k) grid position of objective
    firstOffMax, lastZero, jump2 some counters
*    numk(k) ordinal value of k starting with 1
    numg(g) ordinal value of g starting with 0
    step(k) step of grid points in objective functions
    firstjump(k) jumps in the grid points' traversing  - first jumps when all posg(k)=0 - for all objective functions
    firsttime  parameter to perfromr the action only in the very first run
    jump(k) jumps in the grid points' traversing only for the first objective function
;
lastZero=1; loop(km1, numk(km1)=lastZero; lastZero=lastZero+1); numg(g) = ord(g)-1;

grid(km1,g) = yes; // Here we could define different grid intervals for different objectives
maxg(km1) = smax(grid(km1,g), numg(g));
step(km1)=(maxobj(km1)- minobj(km1))/maxg(km1);
gridrhs(grid(km1,g))$(dir(km1)=-1) = maxobj(km1) - numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1));
gridrhs(grid(km1,g))$(dir(km1)=1) = minobj(km1) + numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1));
*display gridrhs;

PUT fx ' Grid points'/   ;
loop (g,
        loop(km1, put gridrhs(km1,g):12:2);
        put /;
     );
put fx /;
put fx 'Efficient solutions'/;

* Walk the grid points and take shortcuts if the model becomes infeasible
posg(km1) = 0;
iter=0;
infeas=0;
firsttime=1;
firstjump(km1)=1;
start=jnow;

repeat
  rhs(km1) = sum(grid(km1,g)$(numg(g)=posg(km1)), gridrhs(km1,g));
  solve mod_epsmethod maximizing a_objval using mip;
  iter=iter+1;
  if (mod_epsmethod.modelstat<>1 and mod_epsmethod.modelstat<>8,  // not optimal is in this case infeasible
    infeas=infeas+1;
    put fx iter:5:0, '  infeasible'/;
    lastZero = 0; loop(km1$(posg(km1)>0 and lastZero=0), lastZero=numk(km1));
    posg(km1)$(numk(km1)<=lastZero) = maxg(km1); // skip all solves for more demanding values of rhs(km1)
  else
    put fx iter:5:0;
    loop(k, put fx z.l(k):12:2);
    put fx ' *** '; // put /;
*
    loop(j, put fx x.l(j):5:0);
    put fx ' *** '; // put /;
*    
    loop(km1, put fx sl.l(km1):12:2, put fx posg(km1):6:0); //put fx sl.l("2"):12:2; put fx sl.l("3"):12:2 ;

    put fx ' *** '; // put /;
    jump(km1)=1;
* only for the first run with all posg(k)=0
    if (firsttime=1,
        loop(km1,
             firstjump(km1)=1+floor(sl.L(km1)/step(km1));
             put fx firstjump(km1):8:0;
             );
        firsttime=0
        );
* calculate only for the first constrained objective function jump(km1)
    put fx ' * ';
*    loop(km1$(numk(km1)=1), jump(km1)=1+floor(sl.L(km1)/step(km1)));
    jump(km1)$(numk(km1)=1)=1+floor(sl.L(km1)/step(km1));
*    put ' lastzero= ', lastzero:3:0 ;
    loop(km1, put fx jump(km1):5:0) ;
    loop(km1$(jump(km1)> 1),  put '   jump')
    put /;
    );
* Proceed forward in the grid
  firstOffMax = 0;
  loop(km1$(posg(km1)<maxg(km1) and firstOffMax=0), posg(km1)=min((posg(km1)+jump(km1)),maxg(km1)); firstOffMax=numk(km1));
* check if the next objective has firstjump
  jump2=0;
  loop(km1,
      if (firstOffMax>1 and numk(km1)=firstOffMax and posg(km1)<=firstjump(km1)-1,jump2=1)
       );
  loop(km1,
      put fx firstOffmax:6, firstjump(km1):6 ;
      if (jump2=1,
           posg(km1)$(numk(km1)<firstOffMax) = firstjump(km1);
*           put fx numk(km1):6, posg(km1):6, '     jump2=1' /;
      else
           posg(km1)$(numk(km1)<firstOffMax) = 0
*           put fx numk(km1):6, posg(km1):6, '    jump2=0' /;
          );
      );
*until sum(km1$(posg(km1)=maxg(km1)),1)= card(km1) and firstOffMax=0;
until sum(km1$(posg(km1)=maxg(km1)),1)= card(km1) and firstOffMax=0;

finish=jnow;
elapsed_time=(finish-start)*86400;

put /;
put 'Infeasibilities = ', infeas:5:0 /;
put 'Elapsed time: ',elapsed_time:10:2, ' seconds' / ;
*$offtext
putclose fx; // close the point file
**$offtext
