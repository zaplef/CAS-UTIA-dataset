SETS
   t   periods       /t1*t4/
   i   products      /V1*V5/
Alias (i,j);

PARAMETERS

p(t)  EUA_price
r(t)  free_allowances
w(i)  production_capacities
h(i)  CO2_coeff
xi    interest_rate
/0.04/

Q  prohibitive_constant
/10000000/

m(i,t) margin
d(i,t) demand
B(i,j) adjusted_technological_matrix:(E-A)^-1;

FREE VARIABLES
             UF        objective_function_value;

POSITIVE VARIABLE
             pp(t)        allowances_purchased
             pm(t)        allowances_sold
             y(i,t)       final_production
             s(t)         banked_allowances
             x(i,t)       raw_production
             E(t)         amount_of_emissions  ;


BINARY VARIABLES
sigma(t)   dummy_binary_variable ;

$CALL GDXXRW.EXE indata.xlsx par=r rng=r!A1 rdim=1 cdim=0 par=p rng=p!A1 rdim=1 cdim=0 par=w rng=w!A1 rdim=1 cdim=0 par=h rng=h!A1 rdim=1 cdim=0 par=m rng=m!A1 rdim=1 cdim=1 par=d rng=d!A1 rdim=1 cdim=1 par=B rng=B!A1 rdim=1 cdim=1
$GDXIN indata.gdx
$LOAD r,p,w,h,m,d,B
$GDXIN

EQUATION

    obj
    de(i,t) demand_const
    C(i,t)  production_cap_const
    P1(t) auxiliary_const1
    P2(t) auxiliary_const2
    y1(t) final_production_iron_const
    vyr(i,t) leontief_model_balance_const
    Em(t)  CO2_emissions_const
    PR1  transfered_allowances_const1
    PR2  transfered_allowances_const2
    PR3  transfered_allowances_const3
    PR4  transfered_allowances_const4
    SP   anti-speculation;


    obj..         UF=e= sum(t,sum(i,y(i,t)*m(i,t))+p(t)*(pm(t)-(1+xi)*pp(t)));
    de(i,t)..     y(i,t)=l=d(i,t);
    vyr(i,t)..    x(i,t)=e=sum(j,B(i,j)*y(j,t));
    C(i,t)..      x(i,t)=l=w(i);
    P1(t)..       pp(t)=l=Q*sigma(t);
    P2(t)..       pm(t)=l=Q*(1-sigma(t));
    y1(t)..       y("V1",t)=e=0;
    Em(t)..       sum(i,x(i,t)*h(i))=e=E(t);
    PR1..         s("t1")=e=r("t1")-E("t1")+pp("t1")-pm("t1");
    PR2..         s("t2")=e=r("t2")-E("t2")+pp("t2")-pm("t2")+s("t1");
    PR3..         s("t3")=e=r("t3")-E("t3")+pp("t3")-pm("t3")+s("t2");
    PR4..         0=e=r("t4")-E("t4")+pp("t4")-pm("t4")+s("t3");
    SP(t)..       pp(t)=l=r(t);


MODEL Fuzzy /all/;

SOLVE Fuzzy USING MIP maximizing UF;
file results/results.txt/;
