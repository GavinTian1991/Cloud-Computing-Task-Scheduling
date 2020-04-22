package aco.ts;

import java.util.Random;

public class aco {
	
	/**
	 * Number of ants
	 */
	int antNum; 
	
	/**
	 * Number of VMs
	 */
	int vmNum;  
	
	/**
	 * Number of tasks
	 */
	int taskNum; 
	
	/**
	 * Number of iterations
	 */
	int iterNum;  
	
	/**
	 *Array of tasks
	 */
	double[] taskArray;  
	
	/**
	 *Array of VMs
	 */
	double[] vmArray;  
	
	/**
	 *Task processing time on VMs
	 */
	double[][] timeMatrix; 
	
	/**
	 *Matrix of pheromone
	 */
	double[][] pheromoneMatrix;
	
	/**
	 * upper and lower boundary of input tasks
	 */
	int taskLengthMax = 100;
	int taskLengthMin = 10;
	
	/**
	 *  upper and lower boundary of VMs processing speed
	 */
	int vmProcessingSpeedMax = 100;
	int vmProcessingSpeedMin = 10;
	
	
	/**
	 * shortest processing time of each ants in each iterations
	 */
	double[] antBestVmTime;
	
	/**
	 * record array of antBestVmTime at each iterations
	 */
	double[][] antBestVmTimeResult;
	
	/**
	 * Number of ant with shortest processing time at each iteration
	 */
	int bestAntTimeNum;
	
	
	/**
	 *Index critical value at each iteration(< use ACO algorithm; > random select)
	 */
	int[] criticalPointMatrix; 
	
	/** 
	 *Index of max pheromone
	 */
	int[] maxPheromoneMatrix;
	
	/** 
	 *path selected by all ant
	 */
	int[][][] pathMatrix_allAnt;
	
	
	/**
	 * alpha
	 */
	private double alpha = 1.0;
	
	/**
	 * beta
	 */
	private double beta = 2.0;
	
	/**
	 * rou
	 */
	private double rou = 0.5;
	
	/**
	 * Q_value
	 */
	double Q_value = 2.0;
	
	//ACO Init
	public void init(int antCount, int vmCount, int iterCount, int taskCount
			,double alphaCount,double betaCount, double rouCount, double qCount) {
		
		this.antNum = antCount;
		this.vmNum = vmCount;
		this.iterNum = iterCount;
		this.taskNum = taskCount;
		this.alpha = alphaCount;
		this.beta = betaCount;
		this.rou = rouCount;
		this.Q_value = qCount;
		
		antBestVmTime = new double[antNum];
		
		antBestVmTimeResult = new double[iterNum][antNum];
		
		pathMatrix_allAnt = new int[antNum][taskNum][vmNum];
		
		taskArray = initRandomArray(taskNum,taskLengthMax,taskLengthMin);
		vmArray = initRandomArray(vmNum,vmProcessingSpeedMax,vmProcessingSpeedMin);
		
		initTimeMatrix(taskArray,vmArray);   
		initPheromoneMatrix(taskNum,vmNum); 
		
		criticalPointMatrix = new int[taskNum];
		maxPheromoneMatrix = new int[taskNum];	
	}
	
	//main ACO algorithm entrance
	public void acoRun() {
		
		//Number of iterations
		for(int gen = 0 ; gen < iterNum;gen++) {
			//Number of ants
			for(int antIndex = 0 ; antIndex < antNum;antIndex++) {
				
				int[][] pathMatrix_oneAnt = new int[taskNum][vmNum];
				//Number of tasks
				for(int taskIndex = 0; taskIndex < taskNum; taskIndex++) {  
					int vmCount = selectVM(antIndex,taskIndex);
					pathMatrix_oneAnt[taskIndex][vmCount] = 1;
				}	
				pathMatrix_allAnt[antIndex] = pathMatrix_oneAnt;
			}
			
			// calculate best processing time of each ant 
			antBestVmTime = calTime_oneIt();
			
			//record best processing time of each ant 
			for(int antIndex = 0 ; antIndex < antNum; antIndex++) {
				antBestVmTimeResult[gen][antIndex] = antBestVmTime[antIndex];
			}
			
			//updatePheromone
			updatePheromone();
			
			//re-init parameters
			reinit();
			
		}
		
	}
	
	public int selectVM(int antCount, int taskCount) {
	    
		//2 ways to select next VMs
	    if (antCount <= criticalPointMatrix[taskCount]) {
	    	//return maxPheromoneMatrix[taskCount];  
	    	return selectVmPara(antCount,taskCount);
	    }
	    else {
	    	//random select VMs
	    	return randomVm(0,vmNum - 1);
	    }	
	}
	
	public int selectVmPara(int antIndex, int taskIndex) {
		
		double[] pVm = new double[vmNum];
		double sum = 0.0;  //信息素概率总和
		
		// Denominator Parts
		for (int vmIndex = 0; vmIndex < vmNum; vmIndex++) {
			//if (selectedVm[i] == 0) {
				sum += Math.pow(pheromoneMatrix[taskIndex][vmIndex], this.alpha)
						* (Math.pow(1.0 / timeMatrix[taskIndex][vmIndex], this.beta));
			//}
		}
		
		//Numerator parts
		for (int vmIndex = 0; vmIndex < vmNum; vmIndex++) {
			//if (selectedVm[i] == 1) {
			//	p[i] = 0.0;
			//} else {
				pVm[vmIndex] = Math.pow(pheromoneMatrix[taskIndex][vmIndex], this.alpha)
						* (Math.pow(1.0 / timeMatrix[taskIndex][vmIndex], this.beta)) / sum;
			//}
		}
		
		int vmselect = getRandomVM(pVm);
		//System.out.println("vmselect = " + vmselect);
		//selectedVm[vmselect] = 1;      
		
		return vmselect;

	}
	
	
	private int getRandomVM(double[] p) {
		double selectP = new Random(System.currentTimeMillis()).nextDouble();
		double sumSel = 0.0;
		for (int i = 0; i < vmNum; i++) {
			sumSel += p[i];
			if (sumSel > selectP)
				return i;
		}
		return -1;
	}
	
	public double[] calTime_oneIt() {
		
		double[] time_allAnt = new double[antNum];
		
		for(int antIndex = 0 ; antIndex < antNum; antIndex++) {
			
			int[][] antBestVmPath = pathMatrix_allAnt[antIndex];
			
			double maxTime = Double.MIN_VALUE;
			
			for(int vmIndex = 0 ; vmIndex < vmNum; vmIndex++) {
				double time = 0.0;
	            for (int taskIndex=0; taskIndex < taskNum; taskIndex++) {
	                if (antBestVmPath[taskIndex][vmIndex] == 1) {
	                    time += timeMatrix[taskIndex][vmIndex];
	                }
	            }
	            if (time > maxTime) {
	                maxTime = time;
	            }
			}
			
//=================Test==================
//			if(antIndex == 30) {
//				System.out.println(maxTime);
//			}
			
			time_allAnt[antIndex] = maxTime;
		}
		return time_allAnt;
		
	}
	
	private void updatePheromone() {
		
		
		for(int i = 0 ; i < taskNum ; i++) {
			for(int j = 0 ; j < vmNum;j++) {
				pheromoneMatrix[i][j] *= rou;
			}
		}
		
	    double minTime = Double.MAX_VALUE;
	    int minAntIndex = -1;
	    for (int antIndex=0; antIndex < antNum; antIndex++) {
	    	if (antBestVmTime[antIndex] < minTime) {
	            minTime = antBestVmTime[antIndex];
	            minAntIndex = antIndex;
	        }
	    	
	    }
	    
	    for (int taskIndex=0; taskIndex<taskNum; taskIndex++) {
	        for (int vmIndex=0; vmIndex<vmNum; vmIndex++) {
	            if (pathMatrix_allAnt[minAntIndex][taskIndex][vmIndex] == 1) {
	                pheromoneMatrix[taskIndex][vmIndex] *= Q_value;
	            }
	        }
	    }
	    
	    //=====================================================================

	    for(int i = 0 ; i < taskNum;i++) {
	    	maxPheromoneMatrix[i] = 0;
	    	criticalPointMatrix[i] = 0;
	    }
	    
	    for (int taskIndex = 0; taskIndex < taskNum; taskIndex++) {
	        
	        int maxIndex = 0;
	        double maxPheromone = pheromoneMatrix[taskIndex][0];
	        double sumPheromone = pheromoneMatrix[taskIndex][0];
	        boolean isAllSame = true;

	        for (int vmIndex=1; vmIndex < vmNum; vmIndex++) {
	            if (pheromoneMatrix[taskIndex][vmIndex] > maxPheromone) {
	                maxPheromone = pheromoneMatrix[taskIndex][vmIndex];
	                maxIndex = vmIndex;
	            }

	            if (pheromoneMatrix[taskIndex][vmIndex] != pheromoneMatrix[taskIndex][vmIndex-1]){
	                isAllSame = false;
	            }

	            sumPheromone += pheromoneMatrix[taskIndex][vmIndex];
	        }

	        if (isAllSame==true) {
	            maxIndex = randomVm(0, vmNum - 1);
	            maxPheromone = pheromoneMatrix[taskIndex][maxIndex];
	        }

	        maxPheromoneMatrix[taskIndex] = maxIndex;	        
	        
	        criticalPointMatrix[taskIndex] = (int)Math.round(antNum * (maxPheromone/sumPheromone));
	        
	}
	    
	    
}
	
	public void reinit() {
		for(int i = 0 ; i < antNum;i++) {
			antBestVmTime[i] = 0;
		}
		
		for(int i = 0 ; i < antNum;i++) {
			for(int j = 0 ; j < taskNum;j++) {
				for(int n = 0 ; n < vmNum;n++) {
					pathMatrix_allAnt[i][j][n] = 0;
				}
			}
		}
	}
	
	
	
	public double[] initRandomArray(int length, int max,int min) {
		
		double[] arrayTemp = new double[length];
		
		for(int i=0;i < length;i++) {
			arrayTemp[i] = min + (int)(Math.random()*(max - min +1));
		}
		
		return arrayTemp;
	}
	
	public void initTimeMatrix(double[] tasks,double[] vms) {
		
		timeMatrix = new double[tasks.length][vms.length];
		
		for(int i = 0 ; i < tasks.length;i++) {
			for(int j = 0 ; j < vms.length;j++) {
				timeMatrix[i][j] = tasks[i]/vms[j];
			}
		}
	}
	
	public void initPheromoneMatrix(int taskNum,int vmNum) {
		
		pheromoneMatrix = new double[taskNum][vmNum];
		for(int i = 0 ; i < taskNum ; i++) {
			for(int j = 0; j < vmNum;j++) {
				pheromoneMatrix[i][j] = 1;
			}
		}
	}
	
	
	public int randomVm(int min ,int max) {
		return (min + (int)(Math.random()*(max - min +1)));
	}
	
	
	
}
