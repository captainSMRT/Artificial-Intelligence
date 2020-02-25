import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Random;
import java.util.Scanner;

public class Perceptron {
	static double[] attribute1  = new double[134];
	static double[] attribute2  = new double[134];
	static double[] attribute3  = new double[134];
	static double[] attribute4  = new double[134];
	static int[] attributeOutput  = new int[134];
	
	static double[] testAttribute1 = new double[16];
	static double[] testAttribute2 = new double[16];
	static double[] testAttribute3 = new double[16];
	static double[] testAttribute4 = new double[16];
	static int[] testAttributeOutput = new int[16];
	
	static double bias = 1;
	static double biasWeight = 1;
	static double attribute1Weight;
	static double attribute2Weight;
	static double attribute3Weight;
	static double attribute4Weight;
	static double error;
	static int match = 0;
	static int termination;
	int epoch = 0;
	double learningRate;
	
	public Perceptron(double learningRate, int termination){
		this.learningRate = learningRate;
		this.termination = termination;
	}
	
	public void beginTraining(){
		attribute1Weight = findWeightByRandom(attribute1);
		attribute2Weight = findWeightByRandom(attribute2);
		attribute3Weight = findWeightByRandom(attribute3);
		attribute4Weight = findWeightByRandom(attribute4);
		
		while(epoch<termination){
			double net = 0;
			for(int i = 0 ; i < attribute1.length; i++){
				int activation;
				net = calculateNet(bias, biasWeight, attribute1[i], attribute1Weight, attribute2[i], attribute2Weight, attribute3[i], attribute3Weight, attribute4[i], attribute4Weight);
				if(net < 0){
					activation = 1;
				}else if(net < 134){
					activation = 2;
				}else{
					activation = 3;
				}
				error = attributeOutput[i] - activation;
				
				if(activation == attributeOutput[i]){
					match = match + 1;
				}
				biasWeight = updateWeight(biasWeight, error, bias,learningRate);
				attribute1Weight = updateWeight(attribute1Weight, error, attribute1[i],learningRate);
				attribute2Weight = updateWeight(attribute2Weight, error, attribute2[i],learningRate);
				attribute3Weight = updateWeight(attribute3Weight, error, attribute3[i],learningRate);
				attribute4Weight = updateWeight(attribute4Weight, error, attribute4[i],learningRate);
				//System.out.println("[Debug] " + (i+1) +" net: " + net);
			}
			System.out.println("Bias weight: " + biasWeight + " | Attribute1Weight: " + attribute1Weight + " | Attribute2Weight: " + attribute2Weight + " | Attribute3Weight: " + attribute3Weight + " | Attribute4Weight: " + attribute4Weight);
			System.out.println("match: " + match);
			if (match!=attribute1.length){
				match = 0;
			}
			epoch = epoch + 1;
		}
		
		System.out.println("Training complete");
		System.out.println("Epoch: " + epoch);
		System.out.println();
		
	}
	
	public void beginTesting(){
		System.out.println("Testing with the following weights...");
		System.out.println("Bias weight: " + biasWeight + " | Attribute1Weight: " + attribute1Weight + " | Attribute2Weight: " + attribute2Weight + " | Attribute3Weight: " + attribute3Weight + " | Attribute4Weight: " + attribute4Weight);
		int counterCorrect = 0;
		int counterWrong = 0;
		for(int i = 0; i < testAttribute1.length; i++){
			int activation = 0;
			double net = calculateNet(bias, biasWeight, testAttribute1[i], attribute1Weight, testAttribute2[i], attribute2Weight, testAttribute3[i], attribute3Weight, testAttribute4[i], attribute4Weight);
			if(net < 0){
				activation = 1;
			}else if(net < 134){
				activation = 2;
			}else{
				activation = 3;
			}
			int match = activation - testAttributeOutput[i];
			boolean isCorrect = false;
			if(match==0){
				counterCorrect = counterCorrect + 1;
				isCorrect = true;
			}else{
				counterWrong = counterWrong + 1;
			}
			String flower;
			if(activation==1){
				flower = "Iris-setosa";
			}else if(activation==2){
				flower = "Iris-versicolor";
			}else{
				flower = "Iris-virginica";
			}
			System.out.println("Testing with element " + (i+1) +": ");
			System.out.println("Match: " + isCorrect + " | Activation: " + activation + " | Target :" + flower + "\n");
		}
		System.out.println("Testing complete");
		System.out.println("Number of matches: " + counterCorrect);
		System.out.println("Number of mismatches: " + counterWrong);

	}
	
	public static void main(String[] args) throws FileNotFoundException {
		FileReader fr = new FileReader("iris.data");
		System.out.println(dataToArrayTraining(fr,134));
		double learningRate = 0.05;
		Perceptron test1 = new Perceptron(learningRate,150000);
		test1.beginTraining();
		System.out.println(attribute1Weight);
		System.out.println(attribute2Weight);
		System.out.println(attribute3Weight);
		System.out.println(attribute4Weight);
		
		FileReader fr2 = new FileReader("iris-testing.data");
		dataToArrayTesting(fr2,16);
		test1.beginTesting();
		

		FileReader fr1 = new FileReader("test.txt");
		dataToArrayTesting(fr1,150);
		test1.beginTesting();
	}
	
	public static boolean dataToArrayTraining(FileReader fr, int size){
		Scanner sc = new Scanner(fr);
		int counter = 0;
		String[] attribute1String = new String[size];
		String[] attribute2String = new String[size];
		String[] attribute3String = new String[size];
		String[] attribute4String = new String[size];
		String[] attributeOutputString = new String[size];
		
		try{
		while(sc.hasNextLine()){
			String[] temp = sc.nextLine().split(",");
			attribute1String[counter] = temp[0];
			attribute2String[counter] = temp[1];
			attribute3String[counter] = temp[2];
			attribute4String[counter] = temp[3];
			attributeOutputString[counter] = temp[4];
			counter = counter + 1;
		}
		
		for(int i = 0; i < attribute1.length;i++){
			attribute1[i] = Double.parseDouble(attribute1String[i]);
			attribute2[i] = Double.parseDouble(attribute2String[i]);
			attribute3[i] = Double.parseDouble(attribute3String[i]);
			attribute4[i] = Double.parseDouble(attribute4String[i]);
			attributeOutput[i] = Integer.parseInt(attributeOutputString[i]);
		}
		
		return true;
		}catch(Exception e){
			e.printStackTrace();
			return false;
		}
	}
	
	public static boolean dataToArrayTesting(FileReader fr, int size){
		
		testAttribute1 = new double[size];
		testAttribute2 = new double[size];
		testAttribute3 = new double[size];
		testAttribute4 = new double[size];
		testAttributeOutput = new int[size];
		
		Scanner sc = new Scanner(fr);
		int counter = 0;
		String[] attribute1String = new String[size];
		String[] attribute2String = new String[size];
		String[] attribute3String = new String[size];
		String[] attribute4String = new String[size];
		String[] attributeOutputString = new String[size];
		
		try{
		while(sc.hasNextLine()){
			String[] temp = sc.nextLine().split(",");
			attribute1String[counter] = temp[0];
			attribute2String[counter] = temp[1];
			attribute3String[counter] = temp[2];
			attribute4String[counter] = temp[3];
			attributeOutputString[counter] = temp[4];
			counter = counter + 1;
		}
		
		for(int i = 0; i < testAttribute1.length;i++){
			testAttribute1[i] = Double.parseDouble(attribute1String[i]);
			testAttribute2[i] = Double.parseDouble(attribute2String[i]);
			testAttribute3[i] = Double.parseDouble(attribute3String[i]);
			testAttribute4[i] = Double.parseDouble(attribute4String[i]);
			testAttributeOutput[i] = Integer.parseInt(attributeOutputString[i]);
		}
		
		return true;
		}catch(Exception e){
			e.printStackTrace();
			return false;
		}
	}
	
	public static double findWeightByRandom(double[] attribute){
		double largest = attribute[0];
		double smallest = attribute[0];
		
		for(int i = 0 ; i < attribute.length; i++){
			if(largest < attribute[i]){
				largest = attribute[i];
			}
			
			if(smallest > attribute[i]){
				smallest = attribute[i];
			}
		}
		
		Random r = new Random();
		int num = r.nextInt((int) ((largest - smallest) + smallest));
		return num + 0.0;
	}
	
	public static double calculateNet(double bias, double w0, double w, double w1, double x, double w2, double y, double w3, double z, double w4){
		return bias * w0 + w * w1 + x * w2 + y * w3 + z * w4;
	}
	
	public static double updateWeight(double weight, double error, double data, double learningRate){
		double newWeight = weight + learningRate * error * data;
		return newWeight;
	}
}


/*
String[] attribute1 = new String[306];
String[] attribute2 = new String[306];
String[] attribute3 = new String[306];

String[] attributeOutput = new String[306];
try{
	
	FileReader fr = new FileReader("haberman.data");
	Scanner sc = new Scanner(fr);
	
	int counter = 0;
	while(sc.hasNextLine()){
		String[] temp = sc.nextLine().split(",");
		attribute1[counter] = temp[0];
		attribute2[counter] = temp[1];
		attribute3[counter] = temp[2];
		attributeOutput[counter] = temp[3];
		counter = counter + 1;
	}
	
}catch(Exception e){
	e.printStackTrace();
}

for(int i = 0; i < attribute1.length; i++){
	System.out.println(attribute1[i] + "," + attribute2[i] + "," + attribute3[i]  + "," + attributeOutput[i]);
}
 
 */