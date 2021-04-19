/*********************************************
 * OPL 12.10.0.0 Model
 * Author: yangk
 * Creation Date: Jun 16, 2020 at 1:38:09 PM
 *********************************************/
 tuple paxTuple {
 	int i;
 	int j;  
 	float t;
 	float pmin;
 	float pmax;
 	float w;
	float l;
 }
 
 tuple rebTuple {
 	int i;
 	int j;  
 	float t;
	float l;	 
 }
 
 tuple edge{
 	int i;
 	int j; 
 }
 
 {string} class = {"auto","non-auto"};
 string path = ...;
 {int} N = ...;
 {paxTuple} PT = ...; 
 {rebTuple} RT = ...;
 float sigma = ...;
 float mu[class] = ...;
 float capacity[class] = ...;
 {edge} P = {<i,j>|<i,j,t,pmin,pmax,w,l> in PT};
 {edge} R = {<i,j>|<i,j,t,l> in RT};
 float maxPrice[P] = [<i,j>:pmax|<i,j,t,pmin,pmax,w,l> in PT];
 float minPrice[P] = [<i,j>:pmin|<i,j,t,pmin,pmax,w,l> in PT];
 float slope[P] = [<i,j>:w|<i,j,t,pmin,pmax,w,l> in PT];
 float paxTime[P] = [<i,j>:t|<i,j,t,pmin,pmax,w,l> in PT];
float paxll[P] = [<i,j>:l|<i,j,t,pmin,pmax,w,l> in PT];
 float rebTime[R] = [<i,j>:t|<i,j,t,l> in RT];
float rebll[R] = [<i,j>:l|<i,j,t,l> in RT];

 dvar float+  q[P];
 dvar float+ x[P][class];
 dvar float+ y[R][class];
 
 execute{
 cplex.predual = 1;
}
 maximize(sum(e in P)q[e]*(maxPrice[e]-slope[e]*q[e]) - sum(e in P,c in class)x[e][c]*(paxll[e]*sigma +mu[c]*paxTime[e])  
 - sum(e in R,c in class)y[e][c]*(rebll[e]*sigma +mu[c]*rebTime[e]));
 
 subject to{
 	forall(i in N)
 	{
 	  forall(c in class)
 	  phi:	sum(e in P:e.i==i)x[e][c] + sum(e in R:e.i==i)y[e][c] ==  sum(e in P:e.j==i)x[e][c] + sum(e in R:e.j==i)y[e][c];
 	} 
 	
 	forall(c in class)
 	gamma:	sum(e in P)x[e][c]*paxTime[e] + sum(e in R)y[e][c]*rebTime[e] <= capacity[c];
 	
 	forall(e in P)
 	{
 	  pi: sum(c in class)x[e][c] == q[e];
 	  alpha:q[e] <= (maxPrice[e]-minPrice[e])/slope[e];
 	}   
 		  	
 }
 
  execute{
	
 	var file = new IloOplOutputFile(thisOplModel.path);
	file.write("v=")
	file.write(gamma["non-auto"].dual+mu["non-auto"]);
 	file.writeln(";")
 	file.write("obj=")
 	file.write(cplex.getObjValue());
 	file.writeln(";")
 	 
 	file.write("x=[")
	for(var e in P)
	{
		file.write("(")
		file.write(e.i);
		file.write(",");
		file.write(e.j);
		file.write(",")
		file.write(x[e]["auto"]);
		file.write(",")
		file.write(x[e]["non-auto"]);
		file.write(",")
		file.write(q[e]);
		file.write(",")
		file.write(q[e]*slope[e]+minPrice[e]);
		file.write(",")
		file.write(-pi[e].dual);
		file.write(")")
	}	
	file.writeln("];");
	
	
	file.write("y=[")
	for(var e in R)
	{
		file.write("(")
		file.write(e.i);
		file.write(",");
		file.write(e.j);
		file.write(",")
		file.write(y[e]["auto"]);
		file.write(",")
		file.write(y[e]["non-auto"]);
		file.write(")")
	}	
	file.writeln("];");
	
	file.write("phi=[")
	for(var i in N)
	{
		file.write("(")
		file.write(i);
		file.write(",");
		file.write(phi[i]["non-auto"].dual);
		file.write(")")
	}	
	file.writeln("];");
 }