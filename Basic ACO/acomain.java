package aco.ts;


public class acomain {
	
	public static void main(String[] args) {	
		
		aco Aco = new aco();
		//antnum vmnum iternum tasknum alpha beta rou Q_value
		Aco.init(40,10,100,100,1.0,2.0,0.5,2.0); 
		Aco.acoRun();
		
	}

}
