SET bus /1*34/
EVloc(bus) /19, 2, 14, 31,28,9,11,16,23,34/
station(bus) /5,20,24/
;
alias(bus,node);
table position(bus,*)
    x       y
1   0.00    14.40
2   3.00    14.20
3   5.00    13.80
4   8.00    14.40
5   12.00   14.40
6   14.50   15.00
7   5.00    12.70
8   5.00    11.20
9   8.00    11.50
10  9.60    11.20
11  12.00   11.00
12  1.00    10.00
13  3.00    9.80
14  4.80    8.50
15  8.00    9.60
16  11.20   7.00
17  4.00    7.00
18  5.50    7.00
19  0.00    4.80
20  2.00    6.00
21  5.00    5.00
22  7.60    5.00
23  9.60    5.00
24  12.80   4.90
25  14.40   4.90
26  3.00    3.00
27  4.80    3.00
28  7.60    3.00
29  11.50   3.10
30  1.50    1.50
31  6.00    1.50
32  9.00    1.50
33  12.30   1.50
34  14.40   1.80;


Set conex 'bus connectivity matrix'
/
   1.2
   
   1.12
   2.3
   2.7
   3.8
   3.4
   4.9
   4.5
   5.10
   5.6
   6.11
   6.25
   7.8
   7.13
   8.9
   8.14
   9.10
   9.15
   10.5
   10.11
   11.16
   11.24
   12.13
   12.19
   13.14
   13.17
   14.15
   14.18
   15.16
   16.23
   17.18
   17.20
   18.21
   19.20
   19.30
   20.26
   20.21
   21.22
   21.27
   22.23
   22.28
   23.24
   24.25
   24.29
   25.34
   26.27
   27.28
   27.31
   28.29
   28.32
   29.34
   29.33
   30.31
   31.32
   32.33
   33.34 
/;
* -----------------------------------------------------
conex(bus,node)$(conex(node,bus)) = 1;



display position;
parameter d(bus,node),demand(bus);
d(bus,node)= sqrt( power(position(bus,'x')-position(node,'x'),2)+power(position(bus,'y')-position(node,'y'),2)  );
display d;

positive variables flow(bus,node),gen(bus);
variable OF;
flow.up(bus,node)=1;

parameter demand(bus);

equations

const1
const2;

const1 .. OF =e=sum( (bus,node),d(bus,node)*flow(bus,node) );
const2(bus) .. gen(bus)-demand(bus) =e=sum(node$conex(bus,node), flow(bus,node) ) - sum(node$conex(bus,node), flow(node,bus) );


model pathfinder /all/;

demand(bus) =0;


alias(EVloc,Evi);
alias(station,si);

parameter report(EVI,SI);
scalar iteration /1/;

set counter /1*30/;


Gen.up(bus)=0;
Gen.up(bus)$(ord(bus) = 19) = 1;
demand(bus)$(ord(bus) = 5) = 1;
solve pathfinder min OF using LP;
display flow.l;



*$ontext

parameter report2(bus,node,counter);
*parameter report3(counter,*);

loop((EVI,SI),

    demand(bus) =0;
    gen.up(bus) =0;

    demand(SI)  =1;
    gen.up(EVI) =1;

    solve pathfinder min OF using LP;
    report(EVI,SI) = OF.l;
    
    report2(bus,node,counter)$(ord(counter)=iteration and flow.l(bus,node)>0)=1 ;
    iteration = iteration +1;
);

display report,report2; 

*$offtext

*report3(counter,'EV')=ord(EVI)$(ord(counter)=iteration) ;
*report3(counter,'Station')=ord(SI)$(ord(counter)=iteration) ;
    






 

