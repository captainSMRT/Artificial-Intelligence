import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class OngKaiXuan_160268922_CO3311cw2 {
	
	public static void main(String[] args){
		
		//Data Set from Course Work 1
		double[] dataX = { -0.82, -0.77, -0.73, -0.71, -0.69, -0.68, -0.66, -0.61, 0.28, 0.31, 0.35, 0.36, 0.37, 0.43,
				0.43, 0.47, 0.47, 0.48, 0.53, 0.54, 0.58, 0.58, 0.65, 0.74 };
		double[] dataY = { 0.49, 0.47, 0.55, 0.61, 0.59, 0.60, 0.58, 0.69, 0.96, 0.95, -0.06, 0.91, 0.93, 0.83, 0.90,
				-0.10, -0.07, 0.82, -0.23, -0.21, -0.38, 0.81, -0.04, -0.05 };
		double[] dataZ = { 0.29, 0.43, 0.41, 0.36, 0.42, 0.42, 0.48, 0.40, 0.01, 0.07, 0.94, -0.22, 0.09, -0.34, 0.01,
				0.88, 0.88, -0.30, 0.82, 0.82, 0.72, 0.00, 0.76, 0.67 };
	

		//Course Work 2 Testing Set
		double[] testingDataX = {-0.99, -0.98, -0.96, -0.94, -0.84, -0.77, -0.73, -0.70, 0.31, 0.31, 0.39, 0.41, 0.46, 0.47, 0.49, 0.50, 0.52, 0.56, 0.57, 0.68, 0.71, 0.80, 0.83, 0.84};
		double[] testingDataY = {0.13, 0.20, 0.29, 0.29, 0.42, 0.41, 0.57, 0.47, -0.05, -0.15, -0.30, -0.14, 0.89, 0.87, 0.87, -0.21, -0.28, -0.35, -0.22, 0.73, 0.69, 0.58, 0.56, 0.54};
		double[] testingDataZ = {0.08, 0.01, 0.03, 0.18, 0.34, 0.48, 0.37, 0.53, 0.95, 0.94, 0.87, 0.90, 0.04, -0.16, 0.09, 0.84, 0.80, 0.75, 0.79, 0.12, 0.13, 0.16, -0.02, 0.10};
				
		//Sample weights if user want to choose selectDefinedWeight(...) instead of selectRandomWeight()
		//double[] unit1Weights = {-0.34, -0.01, 0.77 };  
		//double[] unit2Weights = {-0.01, 0.21, -0.20}; 
		//double[] unit3Weights = {0.45, 0.54};
		
		
		//2 Cluster Training and Testing
		int clusterSize = 2;
		int[] expectedCluster = {2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,1,1,2,1,1,1,2,1,1}; //Classification obtained from Kohonen network of Course Work 1
		BackPropagationTraining bp = new BackPropagationTraining(dataX, dataY, dataZ, expectedCluster, clusterSize);
		bp.setLearningRate(0.1);
		bp.selectRandomWeight();
		//bp.selectDefinedWeight(unit1Weights, unit2Weights, unit3Weights); 
		bp.trainNetwork();
		bp.testNetwork(testingDataX, testingDataY, testingDataZ);
		bp.writeToFile();
		
		//3 Cluster Training and Testing
		clusterSize = 3;
		expectedCluster = new int[]{1,1,1,1,1,1,1,1,3,3,2,3,3,3,3,2,2,3,2,2,2,3,2,2}; //Classification obtained from Kohonen network of Course Work 1
		BackPropagationTraining bp2 = new BackPropagationTraining(dataX, dataY, dataZ, expectedCluster, clusterSize);
		bp2.setLearningRate(0.1);
		bp2.selectRandomWeight();
		//bp2.selectDefinedWeight(unit1Weights, unit2Weights, unit3Weights);
		bp2.trainNetwork();
		bp2.testNetwork(testingDataX, testingDataY, testingDataZ);
		bp2.writeToFile();
		
		//4 Cluster Training and Testing
		clusterSize = 4;
		expectedCluster = new int[]{4,4,4,4,4,4,4,4,3,3,1,3,3,3,3,1,1,3,1,1,1,3,1,1}; //Classification obtained from Kohonen network of Course Work 1
		BackPropagationTraining bp3 = new BackPropagationTraining(dataX, dataY, dataZ, expectedCluster, clusterSize);
		bp3.setLearningRate(0.1);
		bp3.selectRandomWeight();
		//bp3.selectDefinedWeight(unit1Weights, unit2Weights, unit3Weights); 
		bp3.trainNetwork();
		bp3.testNetwork(testingDataX, testingDataY, testingDataZ);
		bp3.writeToFile();
	}
	
}

class BackPropagationTraining{
	
	private int clusterSize;
	private double learningRate;
	private int epoch = 1;
	private int iteration = 1;
	
	private double bias = 1;
	private double [] x;
	private double [] y;
	private double [] z;
	private int [] expectedCluster; //Classification obtained from Kohonen network of Course Work 1
	
	private int [] predictedCluster; //Classification predicted by Backpropagation network
	private int numberOfMatch = 0;  //Count how many matches exist between expectedCluster and predictedCluster
	
	private double unit1WeightBias = 0.5;
	private double unit1WeightX;
	private double unit1WeightY;
	private double unit1WeightZ;
	private double unit1Net;
	private double unit1Activation;
	private double unit1Error;
	
	private double unit2WeightBias = 0.5;
	private double unit2WeightX;
	private double unit2WeightY;
	private double unit2WeightZ;
	private double unit2Net;
	private double unit2Activation;
	private double unit2Error;
	
	private double unit3WeightBias = 1;
	private double unit3WeightUnit1;
	private double unit3WeightUnit2;
	private double unit3Net;
	private double unit3Activation;
	private double [] target;
	private double unit3Error;
	
	//Variables for file writing
	private ArrayList<String> outputLog; //Save all the important details that are to be logged into fullTrainingAndTestingResult
	private FileWriter fullTrainingAndTestingResult; //File will be saved as: clusterSize+"ClusterBP.txt"
	private FileWriter weightChangeIterationResult;//File will be saved as: clusterSize+"SumOfWeightChangePerIteration.txt"
	private FileWriter weightChangeEpochResult; //File will be saved as: clusterSize+"SumOfWeightChangePerEpoch.txt"
	private FileWriter fullWeightsAfterEachEpochResult;//File will be saved as: clusterSize+"AllEpochWeights.txt"
	private ArrayList<ArrayList<Double>> allGeneratedWeightsForEachEpoch = new ArrayList<ArrayList<Double>>(); //Save results into ArrayList for FileWriter
	private ArrayList<ArrayList<Double>> allSumOfWeightsDetails = new ArrayList<ArrayList<Double>>();  //Save results into ArrayList for FileWriter

	public BackPropagationTraining(double[] x, double [] y, double [] z, int [] expectedCluster, int clusterSize){
		
		//Stop the program if either the length of x, y or z doesn't match
		if(x.length!=z.length || y.length!=z.length || x.length!=expectedCluster.length){
			System.out.println("Program shutting down as length of X, Y, Z and expected cluster does not match...");
			System.exit(0);
		}
		
		this.x = x;
		this.y = y;
		this.z = z;
		this.expectedCluster = expectedCluster;
		this.clusterSize = clusterSize;
		
		//Calculation of the targets for Unit 3
		predictedCluster = new int[expectedCluster.length];
		target = new double[x.length];
		double[] range = new double[this.clusterSize+1];
		range[0] = 0;
		for(int i = 1; i < range.length; i++){
			range[i] =  1.0 * i /this.clusterSize;
		}
		double[] calculatedTarget = new double[this.clusterSize];
		for(int i = 0; i < calculatedTarget.length; i++){
			calculatedTarget[i] = (range[i] + range[i+1])/2;
		}
		
		for(int i = 0; i < expectedCluster.length; i++){
			for(int j = 1; j <= this.clusterSize; j++){
				if(expectedCluster[i]==j){
					target[i] = calculatedTarget[j-1];
				}
			}
		}
		
		outputLog = new ArrayList<String>();
		outputLog.add("Cluster size: " + clusterSize);
	}
	
	public void selectRandomWeight(){
		System.out.println("Random weight selected. Initialising weight...");

		//Finding the upper bound and lower bound for all x,y and z so that the respective weights will be values between them
		Random r = new Random();
		double upperBoundX = x[0];
		double lowerBoundX = x[0];
		double upperBoundY = y[0];
		double lowerBoundY = y[0];
		double upperBoundZ = z[0];
		double lowerBoundZ = z[0];
		for(int i = 1 ; i < x.length; i++){
			if(upperBoundX < x[i]){
				upperBoundX = x[i];
			}
			if(upperBoundY < y[i]){
				upperBoundY = y[i];
			}
			if(upperBoundZ < z[i]){
				upperBoundZ = z[i];
			}
			if(lowerBoundX > x[i]){
				lowerBoundX = x[i];
			}
			if(lowerBoundY > y[i]){
				lowerBoundY = y[i];
			}
			if(lowerBoundZ > z[i]){
				lowerBoundZ = z[i];
			}
		}
		
		//Calculating random weights
		unit1WeightX = ( r.nextInt((int)(upperBoundX*100) - (int)(lowerBoundX*100) + 1) + (int)(lowerBoundX*100) ) /100.0;
		unit1WeightY = ( r.nextInt((int)(upperBoundY*100) - (int)(lowerBoundY*100) + 1) + (int)(lowerBoundY*100) ) /100.0;
		unit1WeightZ = ( r.nextInt((int)(upperBoundZ*100) - (int)(lowerBoundZ*100) + 1) + (int)(lowerBoundZ*100) ) /100.0;
		unit2WeightX = ( r.nextInt((int)(upperBoundX*100) - (int)(lowerBoundX*100) + 1) + (int)(lowerBoundX*100) ) /100.0;
		unit2WeightY = ( r.nextInt((int)(upperBoundY*100) - (int)(lowerBoundY*100) + 1) + (int)(lowerBoundY*100) ) /100.0;
		unit2WeightZ = ( r.nextInt((int)(upperBoundZ*100) - (int)(lowerBoundZ*100) + 1) + (int)(lowerBoundZ*100) ) /100.0;
		unit3WeightUnit1 = ( r.nextInt((int)(upperBoundX*100) - (int)(lowerBoundX*100) + 1) + (int)(lowerBoundX*100) ) /100.0;
		unit3WeightUnit2 = ( r.nextInt((int)(upperBoundY*100) - (int)(lowerBoundY*100) + 1) + (int)(lowerBoundY*100) ) /100.0;

		outputLog.add("Random weight selected");
		outputLog.add("Unit 1 Weight X: " + unit1WeightX);
		outputLog.add("Unit 1 Weight Y: " + unit1WeightY);
		outputLog.add("Unit 1 Weight Z: " + unit1WeightZ);
		outputLog.add("Unit 2 Weight X: " + unit2WeightX);
		outputLog.add("Unit 2 Weight Y: " + unit2WeightY);
		outputLog.add("Unit 2 Weight Z: " + unit2WeightZ);
		outputLog.add("Unit 3 Weight for Unit 1's Activation: " + unit3WeightUnit1);
		outputLog.add("Unit 3 Weight for Unit 2's Activation: " + unit3WeightUnit2);
		System.out.println("Random weight initialised");
	}
	
	public void selectDefinedWeight(double[] unit1Weights, double[] unit2Weights, double[] unit3Weights){
		System.out.println("Self defined weights selected. Initialising weight...");
		unit1WeightX = unit1Weights[0];
		unit1WeightY = unit1Weights[1];
		unit1WeightZ = unit1Weights[2];
		
		unit2WeightX = unit2Weights[0];
		unit2WeightY = unit2Weights[1];
		unit2WeightZ = unit2Weights[2];
		
		unit3WeightUnit1 = unit3Weights[0];
		unit3WeightUnit2 = unit3Weights[1];
		
		outputLog.add("Self defined weight selected");
		outputLog.add("Unit 1 Weight X: " + unit1WeightX);
		outputLog.add("Unit 1 Weight Y: " + unit1WeightY);
		outputLog.add("Unit 1 Weight Z: " + unit1WeightZ);
		outputLog.add("Unit 2 Weight X: " + unit2WeightX);
		outputLog.add("Unit 2 Weight Y: " + unit2WeightY);
		outputLog.add("Unit 2 Weight Z: " + unit2WeightZ);
		outputLog.add("Unit 3 Weight for Unit 1's Activation: " + unit3WeightUnit1);
		outputLog.add("Unit 3 Weight for Unit 2's Activation: " + unit3WeightUnit2);
		System.out.println("Self defined weight initialised");
	}
	
	public void trainNetwork(){
		System.out.println("Training condition: Terminate when all predicted cluster matches expected cluster");
		System.out.println("Training network...");
		
		
		outputLog.add("\nPre-training weights");
		outputLog.add("Unit 1 Weight[bias, x, y, z]: " + unit1WeightBias +", " + unit1WeightX + ", " + unit1WeightY + ", " + unit1WeightZ);
		outputLog.add("Unit 2 Weight[bias, x, y, z]: " + unit2WeightBias +", " + unit2WeightX + ", " + unit2WeightY + ", " + unit2WeightZ);
		outputLog.add("Unit 3 Weight[bias, U1, U2]: " + unit3WeightBias +", " + unit3WeightUnit1 + ", " + unit3WeightUnit2);
		
		//SumOfWeightChange is set to 0 because there is the initial weights, and there are no previous set of weights to calculate the changes
		double sumOfWeights = 0.0;
		ArrayList<Double> details = new ArrayList<Double>();
		details.add(0.0);
		details.add(sumOfWeights);
		allSumOfWeightsDetails.add(details);
		addToEpochDetails();
		
		while(numberOfMatch!=x.length){ //This condition will only end if all expectedCluster are equally matched with predictedCluster
			
			//Hold current weight values in preparation to calculate the sum of weight change
			double tempUnit1WeightBias = unit1WeightBias;
			double tempUnit1WeightX = unit1WeightX;
			double tempUnit1WeightY = unit1WeightY;
			double tempUnit1WeightZ = unit1WeightZ;
			double tempUnit2WeightBias = unit2WeightBias;
			double tempUnit2WeightX = unit2WeightX;
			double tempUnit2WeightY = unit2WeightY;
			double tempUnit2WeightZ = unit2WeightZ;
			double tempUnit3WeightBias = unit3WeightBias;
			double tempUnit3WeightUnit1 = unit3WeightUnit1;
			double tempUnit3WeightUnit2 = unit3WeightUnit2;
			
			outputLog.add("================ Epoch " + epoch + " ================");
			numberOfMatch = 0;
			for(int i = 0 ; i < x.length; i++){
				
				//Calculation of Net and Activation
				unit1Net = calculateNet(bias, unit1WeightBias, x[i], unit1WeightX, y[i], unit1WeightY, z[i], unit1WeightZ);
				unit1Activation = calculateSigmoidActivation(unit1Net);
				unit2Net = calculateNet(bias, unit2WeightBias, x[i], unit2WeightX, y[i], unit2WeightY, z[i], unit2WeightZ);
				unit2Activation = calculateSigmoidActivation(unit2Net);
				unit3Net = calculateNet(bias, unit3WeightBias, unit1Activation, unit3WeightUnit1, unit2Activation, unit3WeightUnit2, 0, 0);
				unit3Activation = calculateSigmoidActivation(unit3Net);
			
				//Calculation of errors
				unit3Error = calculateOutputError(target[i], unit3Activation);
				unit1Error = calculateHiddenUnitError(unit3Error, unit1Activation, unit3WeightUnit1);
				unit2Error = calculateHiddenUnitError(unit3Error, unit2Activation, unit3WeightUnit2);
			
				//Calculation of new sets of weights
				unit3WeightBias =  calculateNewWeight(unit3WeightBias,  unit3Error, bias);
				unit3WeightUnit1 =  calculateNewWeight(unit3WeightUnit1,  unit3Error, unit1Activation);
				unit3WeightUnit2 =  calculateNewWeight(unit3WeightUnit2,  unit3Error, unit2Activation);
				unit1WeightBias = calculateNewWeight(unit1WeightBias,  unit1Error, bias);
				unit1WeightX = calculateNewWeight(unit1WeightX,  unit1Error, x[i]);
				unit1WeightY = calculateNewWeight(unit1WeightY,  unit1Error, y[i]);
				unit1WeightZ = calculateNewWeight(unit1WeightZ,  unit1Error, z[i]);
				unit2WeightBias = calculateNewWeight(unit2WeightBias,  unit2Error, bias);
				unit2WeightX = calculateNewWeight(unit2WeightX,  unit2Error, x[i]);
				unit2WeightY = calculateNewWeight(unit2WeightY,  unit2Error, y[i]);
				unit2WeightZ = calculateNewWeight(unit2WeightZ,  unit2Error, z[i]);
				
				//Calculating the sum of weight change
				sumOfWeights = Math.pow((unit1WeightBias - tempUnit1WeightBias),2)  + Math.pow((unit1WeightX - tempUnit1WeightX),2) + Math.pow((unit1WeightY - tempUnit1WeightY),2) + Math.pow((unit1WeightZ - tempUnit1WeightZ),2) +
						Math.pow((unit2WeightBias - tempUnit2WeightBias),2)  + Math.pow((unit2WeightX - tempUnit2WeightX),2) + Math.pow((unit2WeightY - tempUnit2WeightY),2) + Math.pow((unit2WeightZ - tempUnit2WeightZ),2) +
						Math.pow((unit3WeightBias - tempUnit3WeightBias),2)  + Math.pow((unit3WeightUnit1 - tempUnit3WeightUnit1),2) + Math.pow((unit3WeightUnit2 - tempUnit3WeightUnit2),2);
				
				//Saving iteration results that are to be logged in weightChangeIterationResult txt file
				ArrayList tempDetails = new ArrayList();
				tempDetails.add(iteration);
				tempDetails.add(sumOfWeights);
				allSumOfWeightsDetails.add(tempDetails);
				
				//Saving results that are to be logged in fullTrainingAndTestingResult txt file
				outputLog.add("Iteration: " + iteration);
				outputLog.add("Unit 1 Weight[bias, x, y, z]: " + unit1WeightBias +", " + unit1WeightX + ", " + unit1WeightY + ", " + unit1WeightZ);
				outputLog.add("Unit 2 Weight[bias, x, y, z]: " + unit2WeightBias +", " + unit2WeightX + ", " + unit2WeightY + ", " + unit2WeightZ);
				outputLog.add("Unit 3 Weight[bias, U1, U2]: " + unit3WeightBias +", " + unit3WeightUnit1 + ", " + unit3WeightUnit2);
				outputLog.add("Activation: " + unit3Activation);
				predictedCluster[i] = calculatePrediction(unit3Activation); 
				if(predictedCluster[i]==expectedCluster[i]){
					outputLog.add("Predicted Cluster matches Expected Cluster of: " + expectedCluster[i] + "\n");
					numberOfMatch = numberOfMatch + 1;
				}else{
					outputLog.add("Predicted Cluster does not match Expected Cluster. Predicted Cluster: " + predictedCluster[i] +"    |   Expected Cluster: "+ expectedCluster[i]  + "\n");
				}
				
				iteration = iteration + 1;
			}
			
			//Saving epoch results that are to be logged in weightChangeEpochResult txt file
			addToEpochDetails();
			
			outputLog.add("Number of Matches: " + numberOfMatch  + "\n");
			if(numberOfMatch==x.length){ //Saving final training results that are to be logged in fullTrainingAndTestingResult txt file
				iteration = iteration - 1;
				outputLog.add("Training completed. Epoch: " + epoch +"    |    Iteration: " + iteration);
				System.out.println("Training completed. Epoch: " + epoch +"    |    Iteration: " + iteration);
				System.out.println("Finalized weight: ");
				System.out.println("Unit 1 Weight Bias: " + unit1WeightBias);
				System.out.println("Unit 1 Weight X: " + unit1WeightX);
				System.out.println("Unit 1 Weight Y: " + unit1WeightY);
				System.out.println("Unit 1 Weight Z: " + unit1WeightZ);
				System.out.println("Unit 2 Weight Bias: " + unit2WeightBias);
				System.out.println("Unit 2 Weight X: " + unit2WeightX);
				System.out.println("Unit 2 Weight Y: " + unit2WeightY);
				System.out.println("Unit 2 Weight Z: " + unit2WeightZ);
				System.out.println("Unit 3 Weight Bias: " + unit3WeightBias);
				System.out.println("Unit 3 Weight for Unit 1's Activation: " + unit3WeightUnit1);
				System.out.println("Unit 3 Weight for Unit 2's Activation: " + unit3WeightUnit2);
			}else{
				epoch = epoch + 1;
			}
			
		}
	}
	
	//Method to test the Network 
	public void testNetwork(double[] dataX, double[] dataY, double[] dataZ){
		outputLog.add("\n");
		outputLog.add("Testing network using the following data set: ");
		System.out.println("\n");
		System.out.println("Testing network using the following data set: ");
		
		for(int i = 0 ; i < dataX.length; i++){
			outputLog.add("X: " + dataX[i] +" Y: " + dataY[i] + " Z: " + dataZ[i]);
			System.out.println("X: " + dataX[i] +" Y: " + dataY[i] + " Z: " + dataZ[i]);
		}
		
		outputLog.add("\n");
		outputLog.add("Testing network using the following weights: ");
		System.out.println("\n");
		System.out.println("Testing network using the following weights: ");
		
		outputLog.add("Unit 1 Bias Weight: " + unit1WeightBias + " Unit 1 X Weight: " + unit1WeightX + " Unit 1 Y Weight: " + unit1WeightY + " Unit 1 Z Weight: " + unit1WeightZ);
		outputLog.add("Unit 2 Bias Weight: " + unit2WeightBias + " Unit 2 X Weight: " + unit2WeightX + " Unit 2 Y Weight: " + unit2WeightY + " Unit 2 Z Weight: " + unit2WeightZ);
		outputLog.add("Unit 3 Bias Weight: " + unit3WeightBias + " Unit 3 X Weight: " + unit3WeightUnit1 + " Unit 3 Y Weight: " + unit3WeightUnit2);
		
		System.out.println("Unit 1 Bias Weight: " + unit1WeightBias + " Unit 1 X Weight: " + unit1WeightX + " Unit 1 Y Weight: " + unit1WeightY + " Unit 1 Z Weight: " + unit1WeightZ);
		System.out.println("Unit 2 Bias Weight: " + unit2WeightBias + " Unit 2 X Weight: " + unit2WeightX + " Unit 2 Y Weight: " + unit2WeightY + " Unit 2 Z Weight: " + unit2WeightZ);
		System.out.println("Unit 3 Bias Weight: " + unit3WeightBias + " Unit 3 X Weight: " + unit3WeightUnit1 + " Unit 3 Y Weight: " + unit3WeightUnit2);
		
		
		outputLog.add("\nTesting network...");
		System.out.println("\nTesting network...");
		for(int i = 0; i < dataX.length; i++){
			double netUnit1 = calculateNet(bias, unit1WeightBias, dataX[i], unit1WeightX, dataY[i], unit1WeightY, dataZ[i], unit1WeightZ);
			double activationUnit1 = calculateSigmoidActivation(netUnit1);
			double netUnit2 = calculateNet(bias, unit2WeightBias, dataX[i], unit2WeightX, dataY[i], unit2WeightY, dataZ[i], unit2WeightZ);
			double activationUnit2 = calculateSigmoidActivation(netUnit2);
			double netUnit3 = calculateNet(bias, unit3WeightBias, activationUnit1, unit3WeightUnit1, activationUnit2, unit3WeightUnit2, 0, 0);
			double activationUnit3 = calculateSigmoidActivation(netUnit3);
			
			int predictedCluster = calculatePrediction(activationUnit3); 
			outputLog.add("X: " + dataX[i] +" Y: " + dataY[i] + " Z: " + dataZ[i] +" belongs to cluster " + predictedCluster);
			System.out.println("X: " + dataX[i] +" Y: " + dataY[i] + " Z: " + dataZ[i] +" belongs to cluster " + predictedCluster);
		}
		outputLog.add("\nTesting network completed");
		System.out.println("\nTesting network completed");
	}
	//End of Method to test the Network
	
	public void setLearningRate(double learningRate){
		this.learningRate = learningRate;
		outputLog.add("Learning Rate: " + learningRate);
	}

	public void setClusterSize(int clusterSize){
		this.clusterSize = clusterSize;
		outputLog.add("Cluster Size: " + clusterSize);
	}

	public double calculateNet(double bias, double biasWeight, double x, double xWeight, double y, double yWeight, double z, double zWeight){
		return bias * biasWeight + x * xWeight + y * yWeight + z * zWeight;
	}
	
	public double calculateSigmoidActivation(double net){
		return (1 / ( 1 + Math.exp(-net)));
	}
	
	public double calculateOutputError(double target, double activation){
		return (target-activation) * activation * (1 - activation);
	}
	
	public double calculateHiddenUnitError(double error, double activation, double weight){
		return error * activation * (1 - activation) * weight;
	}
	
	public double calculateNewWeight(double previousWeight, double hiddenUnitError, double previousData){
		return previousWeight + learningRate * hiddenUnitError * previousData;
	}
	
	public int calculatePrediction(double activation){
		for(int i = 0; i < clusterSize; i++){
			double range = (clusterSize-i-1) / (clusterSize * 1.0);
			if(activation>=range){				
				return clusterSize - i;
			}
		}
		System.exit(0); //For safety precaution in case the "if" statement above did not return a value. Goal is to never return -1 but terminates immediately if any anomalies occur here
		return -1;
	}
	
	public void addToEpochDetails(){
		ArrayList<Double> epochDetails = new ArrayList<Double>();
		epochDetails.add(unit1WeightBias);
		epochDetails.add(unit1WeightX);
		epochDetails.add(unit1WeightY);
		epochDetails.add(unit1WeightZ);
		epochDetails.add(unit2WeightBias);
		epochDetails.add(unit2WeightX);
		epochDetails.add(unit2WeightY);
		epochDetails.add(unit2WeightZ);
		epochDetails.add(unit3WeightBias);
		epochDetails.add(unit3WeightUnit1);
		epochDetails.add(unit3WeightUnit2);
		allGeneratedWeightsForEachEpoch.add(epochDetails);
	}
	
	//Method to write to File
	public boolean writeToFile(){
		try {
			fullTrainingAndTestingResult = new FileWriter(clusterSize+"ClusterBP.txt");
			for(int i = 0; i < outputLog.size(); i++){
				fullTrainingAndTestingResult.write(outputLog.get(i) + "\r\n");
			}
			
			weightChangeIterationResult = new FileWriter(clusterSize+"SumOfWeightChangePerIteration.txt");
			weightChangeIterationResult.write("Iteration\tSumOfWeightChangePerIteration\r\n");
			for(int i = 0; i < allSumOfWeightsDetails.size();i++){
				weightChangeIterationResult.write(allSumOfWeightsDetails.get(i).get(0) + "\t" + allSumOfWeightsDetails.get(i).get(1) + "\r\n");
			}
			
			weightChangeEpochResult = new FileWriter(clusterSize+"SumOfWeightChangePerEpoch.txt");
			weightChangeEpochResult.write("Epoch\tSumOfWeightChangePerEpoch\r\n");
			for(int i = 1; i < allGeneratedWeightsForEachEpoch.size(); i++){
				
				//Calculation of the Sum of Weight change between two epoch
				double temp = Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(0) - allGeneratedWeightsForEachEpoch.get(i-1).get(0)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(1) - allGeneratedWeightsForEachEpoch.get(i-1).get(1)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(2) - allGeneratedWeightsForEachEpoch.get(i-1).get(2)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(3) - allGeneratedWeightsForEachEpoch.get(i-1).get(3)),2) + 
						Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(4) - allGeneratedWeightsForEachEpoch.get(i-1).get(4)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(5) - allGeneratedWeightsForEachEpoch.get(i-1).get(5)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(6) - allGeneratedWeightsForEachEpoch.get(i-1).get(6)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(7) - allGeneratedWeightsForEachEpoch.get(i-1).get(7)),2) + 
						Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(8) - allGeneratedWeightsForEachEpoch.get(i-1).get(8)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(9) - allGeneratedWeightsForEachEpoch.get(i-1).get(9)),2) + Math.pow((allGeneratedWeightsForEachEpoch.get(i).get(10) - allGeneratedWeightsForEachEpoch.get(i-1).get(10)),2);
				weightChangeEpochResult.write(i +"\t" + temp+ "\r\n");
			}
			
			fullWeightsAfterEachEpochResult = new FileWriter(clusterSize+"AllEpochWeights.txt");
			fullWeightsAfterEachEpochResult.write("Epoch\tSumOfWeightChangePerEpoch\r\n");
			for(int i = 0; i < epoch; i++){
				fullWeightsAfterEachEpochResult.write("Epoch "+ i +": " + allGeneratedWeightsForEachEpoch.get(i));
				fullWeightsAfterEachEpochResult.write("\r\n");
			}
			
			fullTrainingAndTestingResult.close();	
			weightChangeIterationResult.close();
			weightChangeEpochResult.close();
			fullWeightsAfterEachEpochResult.close();
			
		} catch (IOException e) {
			System.out.println("Write to file fail. Terminating program...");
			e.printStackTrace();
			return false;
		}
		return true;
	}
	//End of Method to write to File
	
	
}


