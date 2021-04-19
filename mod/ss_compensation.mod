/*********************************************
 * OPL 12.10.0.0 Model
 * Author: yangk
 * Creation Date: Jun 16, 2020 at 1:38:09 PM
 *********************************************/
 tuple paxTuple {
 	int i;
 	int j;  
 	float t;
 	float p;
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
 {edge} P = {<i,j>|<i,j,t,p,w,l> in PT};
 {edge} R = {<i,j>|<i,j,t,l> in RT};
 float demand[P] = [<i,j>:w|<i,j,t,p,w,l> in PT];
 float price[P] = [<i,j>:p|<i,j,t,p,w,l> in PT];
 float paxTime[P] = [<i,j>:t|<i,j,t,p,w,l> in PT];
 float paxll[P] = [<i,j>:l|<i,j,t,p,w,l> in PT];
 float rebTime[R] = [<i,j>:t|<i,j,t,l> in RT];
 float rebll[R] = [<i,j>:l|<i,j,t,l> in RT];

 dvar float+ x[P][class];
 dvar float+ y[R][class];
 maximize(sum(e in P,c in class)x[e][c]*price[e] - sum(e in P,c in class)x[e][c]*(paxll[e]*sigma +mu[c]*paxTime[e]) 
 - sum(e in R,c in class)y[e][c]*(rebll[e]*sigma +mu[c]*rebTime[e])) ;
 
 subject to{
 	forall(i in N)
 	{
 	  forall(c in class)
 		phi:sum(e in P:e.i==i)x[e][c] + sum(e in R:e.i==i)y[e][c] ==  sum(e in P:e.j==i)x[e][c] + sum(e in R:e.j==i)y[e][c];
 	} 
 	
 	forall(c in class)
 	gamma: sum(e in P)x[e][c]*paxTime[e] + sum(e in R)y[e][c]*rebTime[e] <= capacity[c];
 	
 	forall(e in P)
 	{
 	pi:  sum(c in class)x[e][c] <= demand[e];
 	}   
 		  	
 }
 
  execute{
 	var file = new IloOplOutputFile(thisOplModel.path);
	file.write("v=")
	file.write(gamma["non-auto"].dual+mu["non-auto"]);
 	file.writeln(";")
	file.write("gamma=")
	file.write(gamma["non-auto"].dual);
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
		file.write(x[e]["auto"]+x[e]["non-auto"]);
		file.write(",")
		file.write((price[e]-pi[e].dual));
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
		file.write(-phi[i]["non-auto"].dual);
		file.write(")")
	}	
	file.writeln("];");
 }