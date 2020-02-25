import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

public class OngKaiXuan_160268922_CO3311cw1 {
    public static void main(String[] args) {
/*		//Testing with AND gate        
        double[] bias = {1,1,1,1};
        double[] x = {0, 1, 0, 1};
        double[] y = {0, 0, 1, 1};

        double[] unitBias = {1, -2};
        double[] unitX = {-1, 1};
        double[] unitY = {-1, 1};
        int numberOfUnit = 2;
        double zeroLearningRate = 0;
        int oneEpoch = 1;
        KohonenNetwork networkA = new KohonenNetwork(bias, x, y, numberOfUnit, zeroLearningRate, oneEpoch);
        networkA.selectUnitOfChoice(unitBias, unitX, unitY);
        networkA.trainNetwork();
        networkA.displayTrainingResults();
        networkA.calculateCluster();
        try {
            networkA.fw.close();
            System.out.println();
            System.out.println("output.txt is created. It consist of choice of unit selection, data set, unit, normalised data set, normalised unit and all training results");
            System.out.println("trainingresult.txt is created. It consist of post-training, pre-training and summary of weight change for each epochs");
            System.out.println("cluster.txt is created. It consists of which data in the data set belongs to which units");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
*/    

///*    //Actual training using given data set
    	double[] dataX = { -0.82, -0.77, -0.73, -0.71, -0.69, -0.68, -0.66, -0.61, 0.28, 0.31, 0.35, 0.36, 0.37, 0.43,
                    0.43, 0.47, 0.47, 0.48, 0.53, 0.54, 0.58, 0.58, 0.65, 0.74 };
    	double[] dataY = { 0.49, 0.47, 0.55, 0.61, 0.59, 0.60, 0.58, 0.69, 0.96, 0.95, -0.06, 0.91, 0.93, 0.83, 0.90,
                    -0.10, -0.07, 0.82, -0.23, -0.21, -0.38, 0.81, -0.04, -0.05 };
    	double[] dataZ = { 0.29, 0.43, 0.41, 0.36, 0.42, 0.42, 0.48, 0.40, 0.01, 0.07, 0.94, -0.22, 0.09, -0.34, 0.01,
                    0.88, 0.88, -0.30, 0.82, 0.82, 0.72, 0.00, 0.76, 0.67 };

    	validify(dataX, dataY, dataZ);
        
        int numberOfUnits = 6;
        double learningRate = 0.15;
        int epoch = 60;
        KohonenNetwork network = new KohonenNetwork(dataX, dataY, dataZ, numberOfUnits, learningRate, epoch);

        network.selectRandomUnit();
        network.trainNetwork();
        network.displayTrainingResults();

        network.calculateCluster();
        try {
            network.fw.close();
            System.out.println();
            System.out.println("output.txt is created. It consist of choice of unit selection, data set, unit, normalised data set, normalised unit and all training results");
            System.out.println("trainingresult.txt is created. It consist of post-training, pre-training and summary of weight change for each epochs");
            System.out.println("cluster.txt is created. It consists of which data in the data set belongs to which units");
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }
    

    public static boolean validify(double[] dataX, double[] dataY, double[] dataZ) {
        if (dataX.length == dataY.length && dataX.length == dataZ.length) {
            return true;
        } else {
        	System.out.println("Size of Data does not match.");
    		System.out.println("Program ended");
    		System.exit(0);
            return false;
        }
    }
}


//*/

class KohonenNetwork {
    public double[] dataX;
    public double[] dataY;
    public double[] dataZ;
    public double[] normalizedDataX;
    public double[] normalizedDataY;
    public double[] normalizedDataZ;

    public double[] unitX;
    public double[] unitY;
    public double[] unitZ;
    public double[] normalizedUnitX;
    public double[] normalizedUnitY;
    public double[] normalizedUnitZ;

    public ArrayList<ArrayList<Double>> allNormalizedUnitX;
    public ArrayList<ArrayList<Double>> allNormalizedUnitY;
    public ArrayList<ArrayList<Double>> allNormalizedUnitZ;
    public double learningRate;
    public int epoch;

    public int numberOfUnits;

    public FileWriter fw;
    public KohonenNetwork(double[] dataX, double[] dataY, double[] dataZ, int numberOfUnits, double learningRate, int epoch) {
        this.dataX = dataX;
        this.dataY = dataY;
        this.dataZ = dataZ;
        this.numberOfUnits = numberOfUnits;

        this.learningRate = learningRate;
        this.epoch = epoch;

        this.allNormalizedUnitX = new ArrayList<ArrayList<Double>>();
        this.allNormalizedUnitY = new ArrayList<ArrayList<Double>>();
        this.allNormalizedUnitZ = new ArrayList<ArrayList<Double>>();

        try {
            fw = new FileWriter("output.txt");
            fw.write("Data Set:" + "\n");
            fw.write("(X, Y, Z)" + "\n");
            System.out.println("Data Set:");
            System.out.println("(X,Y,Z)");
            for(int i = 0; i < dataX.length; i++){
                fw.write(dataX[i] + ", " + dataY[i] + ", " + dataZ[i] + "\n");
                System.out.println(dataX[i] + ", " + dataY[i] + ", " + dataZ[i]);
            }
            fw.write("\n");
            System.out.println("\n");
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }


    public void selectUnitOfChoice(double[] unitX, double[]unitY, double[]unitZ){
        if(unitX.length != this.numberOfUnits){
        	System.out.println("Given units size does not match numberOfUnits.");
        	System.out.println("Program ended");
        	System.exit(0);
        }else{
            for(int i = 0 ; i < unitX.length; i++){
                this.unitX = unitX;
                this.unitY = unitY;
                this.unitZ = unitZ;
            }

            try{
                fw.write("Choice of units selection: Unit of Choice" + "\n");
                fw.write("Generated Units:" + "\n");
                fw.write("(X, Y, Z)" + "\n");
                System.out.println("Choice of units selection: Unit of Choice");
                System.out.println("Generated Units:");
                System.out.println("(X, Y, Z)");
                for(int i = 0; i < unitX.length; i++){
                    fw.write(unitX[i] + ", " + unitY[i] + ", " + unitZ[i] + "\n");
                    System.out.println(unitX[i] + ", " + unitY[i] + ", " + unitZ[i]);
                }
                fw.write("\n");
            }catch(Exception e){
                e.printStackTrace();
            }
        }
    }

    public void selectRandomUnit() {
        unitX = initializeRandomUnit(dataX);
        unitY = initializeRandomUnit(dataY);
        unitZ = initializeRandomUnit(dataZ);

        try{
            fw.write("Choice of units selection: Random" + "\n");
            fw.write("Generated Units:" + "\n");
            fw.write("(X, Y, Z)" + "\n");
            System.out.println("Choice of units selection: Random");
            System.out.println("Generated Units:");
            System.out.println("(X, Y, Z)");
            for(int i = 0; i < unitX.length; i++){
                fw.write(unitX[i] + ", " + unitY[i] + ", " + unitZ[i] + "\n");
                System.out.println(unitX[i] + ", " + unitY[i] + ", " + unitZ[i]);
            }
            fw.write("\n");
        }catch(Exception e){
            e.printStackTrace();
        }
    }


    public double[] initializeRandomUnit(double[] dataSet) {


        double[] unitSet = new double[this.numberOfUnits];

        ArrayList<String> randomInitializer = new ArrayList<String>();
        for(int i = 0 ; i < dataSet.length; i++){
            randomInitializer.add(i+"");
        }


        int[] randomNumbers = new int[unitSet.length];
        Random r = new Random();
        int index = 0;
        while(index < randomNumbers.length){
            int rand =  r.nextInt(randomInitializer.size());
            if(randomInitializer.get(rand).equals("N.A.")){

            }else{
                randomNumbers[index] = rand;
                randomInitializer.set(rand,"N.A.");
                index++;
            }
        }

        for(int i = 0 ; i < unitSet.length; i++){
            unitSet[i] = dataSet[randomNumbers[i]];
        }
        return unitSet;
    }

    public void dataNormalization() {
        normalizedDataX = new double[dataX.length];
        normalizedDataY = new double[dataY.length];
        normalizedDataZ = new double[dataZ.length];

        for (int i = 0; i < dataX.length; i++) {
            double length = Math.sqrt((dataX[i] * dataX[i]) + (dataY[i] * dataY[i]) + (dataZ[i] * dataZ[i]));
            normalizedDataX[i] = dataX[i] / length;
            normalizedDataY[i] = dataY[i] / length;
            normalizedDataZ[i] = dataZ[i] / length;
        }
    }

    public void classNormalization() {
        normalizedUnitX = new double[unitX.length];
        normalizedUnitY = new double[unitY.length];
        normalizedUnitZ = new double[unitZ.length];

        ArrayList<Double> tempListX = new ArrayList<Double>();
        ArrayList<Double> tempListY = new ArrayList<Double>();
        ArrayList<Double> tempListZ = new ArrayList<Double>();
        for (int i = 0; i < unitX.length; i++) {
            double length = Math.sqrt((unitX[i] * unitX[i]) + (unitY[i] * unitY[i]) + (unitZ[i] * unitZ[i]));
            normalizedUnitX[i] = unitX[i] / length;
            normalizedUnitY[i] = unitY[i] / length;
            normalizedUnitZ[i] = unitZ[i] / length;
            tempListX.add(normalizedUnitX[i]);
            tempListY.add(normalizedUnitY[i]);
            tempListZ.add(normalizedUnitZ[i]);
        }

        allNormalizedUnitX.add(tempListX);
        allNormalizedUnitY.add(tempListY);
        allNormalizedUnitZ.add(tempListZ);
    }

    public void updateUnit(int epochIndex) {
        int highestNetUnit = 0;
        double highestNet = calculateNet(0, epochIndex);
        for (int i = 0; i < normalizedUnitX.length; i++) {
            double currentNet = calculateNet(i, epochIndex);
            if (currentNet > highestNet) {
                highestNetUnit = i;
                highestNet = currentNet;
            }
        }
        try {
            fw.write("Training using data " + (epochIndex+1) + "[X:" + dataX[epochIndex] + " Y:" + dataY[epochIndex] + " Z:" + dataZ[epochIndex] +"]: " +"\n");
            fw.write("Unit " + (highestNetUnit+1)+" has the highest net." + "\n");
            fw.write("Previous Unit " + (highestNetUnit+1)+": " +"\n");
            fw.write("Previous Unit " + (highestNetUnit+1)+": [" + normalizedUnitX[highestNetUnit] + "  |  " +
                    normalizedUnitY[highestNetUnit] + "  |  " +
                    normalizedUnitZ[highestNetUnit] + "]" + "\n");
            fw.write("Updating unit.." + "\n");
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        System.out.println("Training using data " + (epochIndex+1) + "[X:" + dataX[epochIndex] + " Y:" + dataY[epochIndex] + " Z:" + dataZ[epochIndex] +"]: ");
        System.out.println("Unit " + (highestNetUnit+1)+" has the highest net.");
        System.out.println("Previous Unit " + (highestNetUnit+1)+": [" + normalizedUnitX[highestNetUnit] + "  |  " +
                normalizedUnitY[highestNetUnit] + "  |  " +
                normalizedUnitZ[highestNetUnit] + "]");
        System.out.println("Updating unit..");
        double tempUnitX = normalizedUnitX[highestNetUnit]
                + learningRate * (normalizedDataX[epochIndex] - normalizedUnitX[highestNetUnit]);
        double tempUnitY = normalizedUnitY[highestNetUnit]
                + learningRate * (normalizedDataY[epochIndex] - normalizedUnitY[highestNetUnit]);
        double tempUnitZ = normalizedUnitZ[highestNetUnit]
                + learningRate * (normalizedDataZ[epochIndex] - normalizedUnitZ[highestNetUnit]);
        double length = calculateLength(tempUnitX, tempUnitY, tempUnitZ);

        normalizedUnitX[highestNetUnit] = tempUnitX / length;
        normalizedUnitY[highestNetUnit] = tempUnitY / length;
        normalizedUnitZ[highestNetUnit] = tempUnitZ / length;
        displayNewUnits();


    }

    public double calculateLength(double x, double y, double z) {
        return Math.sqrt(x * x + y * y + z * z);
    }

    public double calculateNet(int unitIndex, int dataIndex) {
        double net = normalizedUnitX[unitIndex] * normalizedDataX[dataIndex]
                + normalizedUnitY[unitIndex] * normalizedDataY[dataIndex]
                + normalizedUnitZ[unitIndex] * normalizedDataZ[dataIndex];
        return net;
    }

    public void trainNetwork() {
        dataNormalization();
        classNormalization();
        try{
            fw.write("Normalized Data Sets: " + "\n");
            fw.write("(X, Y, Z)" + "\n");
            System.out.println("Normalized Data Sets: ");
            System.out.println("(X, Y, Z)");
            for(int i = 0; i < normalizedDataX.length; i++){
                fw.write(normalizedDataX[i] +", " + normalizedDataY[i] + ", " + normalizedDataZ[i] + "\n");
                System.out.println(normalizedDataX[i] +", " + normalizedDataY[i] + ", " + normalizedDataZ[i]);
            }
            System.out.println();
            fw.write("\n");

            fw.write("Normalized Units: " + "\n");
            fw.write("(X, Y, Z)" + "\n");
            System.out.println("Normalized Units: ");
            System.out.println("(X, Y, Z)");
            displayCurrentUnits();

            for (int i = 0; i < epoch; i++) {
                fw.write("################Start of Epoch "+ (i+1) + " ################" + "\n");
                fw.write("Pre-training Normalised Unit Set: " + "\n");
                System.out.println("################Start of Epoch "+ (i+1) + " ################");
                System.out.println("Pre-training Normalised Unit Set: ");

                displayCurrentUnits();
                for (int j = 0; j < normalizedDataX.length; j++) {
                    updateUnit(j);
                }

                fw.write("################End of Epoch "+ (i+1) + " ################" + "\n");
                System.out.println("################End of Epoch "+ (i+1) + " ################");

                ArrayList<Double> tempListX = new ArrayList<Double>();
                ArrayList<Double> tempListY = new ArrayList<Double>();
                ArrayList<Double> tempListZ = new ArrayList<Double>();
                for(int j = 0; j < normalizedUnitX.length; j++){
                    tempListX.add(normalizedUnitX[j]);
                    tempListY.add(normalizedUnitY[j]);
                    tempListZ.add(normalizedUnitZ[j]);
                }
                allNormalizedUnitX.add(tempListX);
                allNormalizedUnitY.add(tempListY);
                allNormalizedUnitZ.add(tempListZ);

            }
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    public void calculateCluster(){
        ArrayList<String> clusters = new ArrayList();
        for(int i = 0; i < this.dataX.length; i++){
            clusters.add("Data " + (i+1) +": " + dataX[i] +", " + dataY[i] + ", "+ dataZ[i]);
        }

        for(int i = 0; i < this.dataX.length; i++){
            double highestNet = calculateNet(0, i);
            int highestUnit = 0;
            for(int j = 0; j < this.numberOfUnits ; j++){
                double temp = calculateNet(j, i);
                if (temp > highestNet){
                    highestNet = temp;
                    highestUnit = j;
                }
            }
            clusters.set(i,clusters.get(i)+" Belongs to Unit " + (highestUnit+1));
        }
        try {
            FileWriter fw = new FileWriter("cluster.txt");
            for(int i = 0 ; i < clusters.size(); i++){
                fw.write(clusters.get(i) +"\n");
                System.out.println(clusters.get(i));
            }
            fw.close();
        } catch (IOException e) {

            e.printStackTrace();
        }


    }
    //All of the below methods are display methods that displays key factors such as Data, Class and training results
    public void displayData() {
        DecimalFormat df = new DecimalFormat("#.############");
        for (int i = 0; i < dataX.length; i++) {
            System.out.println(
                    "[" + df.format(dataX[i]) + "  |  " + df.format(dataY[i]) + "  |  " + df.format(dataZ[i]) + "]");
        }
    }

    public void displayClass() {
        for (int i = 0; i < unitX.length; i++) {
            System.out.println("[" + unitX[i] + "  |  " + unitY[i] + "  |  " + unitZ[i] + "]");
        }
    }


    public void displayCurrentUnits() {
        try{
            for (int i = 0; i < normalizedUnitX.length; i++) {
                fw.write(
                        "Unit " + (i+1) +": [" + normalizedUnitX[i] + "  |  " + normalizedUnitY[i] + "  |  " + normalizedUnitZ[i] + "]"+"\n");
                System.out.println(
                        "Unit " + (i+1) +": [" + normalizedUnitX[i] + "  |  " + normalizedUnitY[i] + "  |  " + normalizedUnitZ[i] + "]");
            }
            fw.write("\n");
            System.out.println();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void displayNewUnits() {
        try{
            for (int i = 0; i < normalizedUnitX.length; i++) {
                fw.write(
                        "Updated Unit " + (i+1) +": [" + normalizedUnitX[i] + "  |  " + normalizedUnitY[i] + "  |  " + normalizedUnitZ[i] + "]"+"\n");
                System.out.println(
                        "Updated Unit " + (i+1) +": [" + normalizedUnitX[i] + "  |  " + normalizedUnitY[i] + "  |  " + normalizedUnitZ[i] + "]");
            }
            fw.write("\n");
            System.out.println();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void displayTrainingResults(){
        try {
            FileWriter trainingResult = new FileWriter("trainingresult.txt");
            trainingResult.write("################################Training Result################################" + "\n");
            for(int i = 0; i < epoch; i++){
                trainingResult.write(">>>>>>>>>>Epoch " + (i+1) + ": " +"\n");

                double overallWeightChange = 0;
                for(int j = 0; j < allNormalizedUnitX.get(i).size(); j++){
                    trainingResult.write("Post-training Unit " + (j+1) + ":\t" + allNormalizedUnitX.get(i+1).get(j) +"\t"
                            + allNormalizedUnitY.get(i+1).get(j) + "\t"
                            + allNormalizedUnitZ.get(i+1).get(j) + "\t" +"\n");
                }

                trainingResult.write("\n");

                for(int j = 0; j < allNormalizedUnitX.get(i).size(); j++){
                    trainingResult.write("Pre-training Unit " + (j+1) + ":\t" + allNormalizedUnitX.get(i).get(j) +"\t"
                            + allNormalizedUnitY.get(i).get(j) + "\t"
                            + allNormalizedUnitZ.get(i).get(j) + "\t" +"\n");
                }

                trainingResult.write("\n");

                for(int j = 0; j < allNormalizedUnitX.get(i).size(); j++){
                    trainingResult.write("Changes in Unit " + (j+1) + ":\t" + (allNormalizedUnitX.get(i+1).get(j) - allNormalizedUnitX.get(i).get(j)) +"\t"
                            + (allNormalizedUnitY.get(i+1).get(j) - allNormalizedUnitY.get(i).get(j)) + "\t"
                            + (allNormalizedUnitZ.get(i+1).get(j) - allNormalizedUnitZ.get(i).get(j)) + "\t" +"\n");

                    overallWeightChange = overallWeightChange + (allNormalizedUnitX.get(i+1).get(j) - allNormalizedUnitX.get(i).get(j)) + (allNormalizedUnitY.get(i+1).get(j) - allNormalizedUnitY.get(i).get(j)) +
                            (allNormalizedUnitZ.get(i+1).get(j) - allNormalizedUnitZ.get(i).get(j));
                }
                trainingResult.write("Overall change of weight:\t" + overallWeightChange+"\n");
                trainingResult.write("\n");
            }
            trainingResult.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
    //End of all display methods
}

