package org.example;

import java.io.*;
import java.util.*;
import java.io.FileWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Main {
    public static final int MAX_ITERATIONS=100;//kmean聚类的最大次数
    public static final int K=3;//每辆车分为3类

    private static final Random r = new Random();
    public static int fe=0;
    public static  int FE=100000;
    public static  int vehicleCount;
    public static int customerCount;
    public static int capacity=300;
    public static int initialSolutionCount=10;//遗传算法群体的大小，也等于蚁群大小
    public static double rho=0.1;//局部蒸发率
    public static double beta=2;//决定了信息素和距离的相对重要性
    public static double q0=0.9;//探索与开发的比例
    public static double alpha=0.1;//全局信息更新
    public static double deltaTau;
    public static double[] bestSolution;
    public static double bestValue;
    public static int lengthIndex;

    public static double length =0;

    public static final int MaxRunTime=10;
    public static int count=0;
    public static void main(String[] args) {
        /*String dataFilePath = "C:\\Users\\86136\\Desktop\\课题\\毕设\\数据\\modify-points\\modify-A-n32-k5.txt"; // 数据文件路径*/
        String locationFile="C:\\Users\\86136\\Desktop\\课题\\毕设\\数据\\realLocations\\n100-c3-ne10-locations.csv";
        String distanceFile="C:\\Users\\86136\\Desktop\\课题\\毕设\\数据\\realLocations\\n100-c3-ne10-distances.csv";
        String saveFilePath = "C:\\Users\\86136\\Desktop\\课题\\毕设\\实验结果\\解集\\n100-c3-ne10.txt";
        String countAndLengthPath="C:\\Users\\86136\\Desktop\\课题\\毕设\\实验结果\\解的数量\\n100-c3-ne10.txt";
        double[][]points=readPoints(locationFile);
        double[][]distanceMatrix=readDistance(distanceFile);
        /*double[][] points = initializeData(dataFilePath);
        double[][] distanceMatrix = calculateDistances(points);*/
        double[] demandArray = extractDemands(points);
       /* vehicleCount = Integer.parseInt(dataFilePath.split("-k")[1].replace(".txt", ""));*/
        vehicleCount=5;
        //vehicleCount=extractVehicleCount(locationFile);
        customerCount=points.length-1;
        int runtime=0;

        while(runtime<MaxRunTime){
            ArrayList<double[]> solutions=initialSolutions(distanceMatrix,demandArray);
            lengthIndex=solutions.get(0).length-1;
            bestSolution=solutions.get(0);
            bestValue=bestSolution[bestSolution.length-1];
            deltaTau=1.0/(bestValue*(customerCount+1));
            double[][] pheromeMatrix=initialPherome(deltaTau);
            sortSolutionsByTotalLength(solutions);
            List<double[]> archive=new ArrayList<>();
            archive=updateArchive(archive,solutions);
            while(fe<FE){
                genetic(solutions,demandArray,distanceMatrix,points);
                ACSSort(solutions,pheromeMatrix,distanceMatrix,archive);
                //solutions=localSearch2(solutions,distanceMatrix,demandArray);
                //localSearch(solutions,distanceMatrix,demandArray);
                sortSolutionsByTotalLength(solutions);
                archive=updateArchive(archive,solutions);

            }
            count+=archive.size();
            length +=bestValue;
            fe=0;
            runtime++;
            saveArchive(archive,saveFilePath);
            saveCountAndLength(countAndLengthPath,archive.size(),archive.get(0)[lengthIndex]);
            // printArchive(archive);
            //System.out.println("satisfy?: "+checkAllSolutionsCapacity(archive,demandArray));
            System.out.println("---------------run"+runtime+"----------------");
            System.out.println("count: "+archive.size());
            System.out.println("shortest length: "+archive.get(0)[lengthIndex]);

        }
        count/=(double)MaxRunTime;
        length/=(double) MaxRunTime;
        System.out.println("--------------"+MaxRunTime+"次平均-----------");
        System.out.println("最短路径长度： "+ length );
        System.out.println("解的数量： "+count);

        saveCountAndLength(countAndLengthPath,count,length);

    }
    public static double[][] readPoints(String locationFile){
        List<double[]> points = new ArrayList<>();
        try(BufferedReader br=new BufferedReader(new FileReader(locationFile))){
             String line;
             br.readLine();
             while((line=br.readLine())!=null){
                 String[] values=line.split(",");
                 double latitude=Double.parseDouble(values[1]);
                 double longtitude=Double.parseDouble(values[2]);
                 double demand=Double.parseDouble(values[3]);
                 points.add(new double[]{latitude,longtitude,demand});
             }
            return points.toArray(new double[0][]);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }
    public static double[][] readDistance(String distanceFile){
        List<double[]> distance=new ArrayList<>();
        try(BufferedReader br=new BufferedReader(new FileReader(distanceFile))){
            String line;
            while((line=br.readLine())!=null){
                String[] values=line.split(",");
                double[] row=Arrays.stream(values).mapToDouble(s->Double.parseDouble(s)).toArray();
                distance.add(row);
            }
            return distance.toArray(new double[0][]);
        }catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    public static int extractVehicleCount(String filePath) {
            // 定义匹配模式：k后面跟一个或多个数字
            Pattern pattern = Pattern.compile("k(\\d+)");
            Matcher matcher = pattern.matcher(filePath);

            // 查找匹配的内容
            if (matcher.find()) {
                // 提取第一个捕获组 (即数字部分)
                return Integer.parseInt(matcher.group(1));
            } else {
                throw new IllegalArgumentException("File path does not contain vehicle count in the expected format.");
            }
    }
    public static void saveCountAndLength(String countAndLengthPath,double count,double length){
        try(BufferedWriter bw=new BufferedWriter(new FileWriter(countAndLengthPath,true))){
            bw.write("count:"+count+",length:"+String.format("%.3f", length));
            bw.newLine();
            bw.close();
        }catch (Exception e){
            System.out.println("写文件出错");
        }
    }
    //初始化数据
    public static double[][] initializeData(String filePath) {
        ArrayList<double[]> pointsList = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;

            while ((line = br.readLine()) != null) {
                String[] values = line.split(" ");
                double[] point = new double[3];
                point[0] = Double.parseDouble(values[1]); // x坐标
                point[1] = Double.parseDouble(values[2]); // y坐标
                point[2] = Double.parseDouble(values[3]); // 需求量
                pointsList.add(point);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // 转换为二维数组
        double[][] points = new double[pointsList.size()][3];
        for (int i = 0; i < pointsList.size(); i++) {
            points[i] = pointsList.get(i);
        }

        return points;
    }
    // 计算欧几里得距离矩阵
    public static double[][] calculateDistances(double[][] points) {
        int n = points.length;
        double[][] distanceMatrix = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double dx = points[i][0] - points[j][0];
                    double dy = points[i][1] - points[j][1];
                    distanceMatrix[i][j] = Math.sqrt(dx * dx + dy * dy); // 计算距离
                } else {
                    distanceMatrix[i][j] = 0; // 同一点距离为0
                }
            }
        }

        return distanceMatrix;
    }
    // 提取需求量数组
    public static double[] extractDemands(double[][] points) {
        int n = points.length;
        double[] demands = new double[n];

        for (int i = 0; i < n; i++) {
            demands[i] = points[i][2]; // 需求量存储在第三列
        }

        return demands;
    }

    //贪心算法得到一个初始解
    public static double[] greedyAlgorithm(double[][] distanceMatrix, double[] demandArray) {
        double[] solution = new double[vehicleCount + customerCount + 2];
        solution[0] = 0;
        int index = 1;
        double totalDistance = 0; // 初始化总路径长度
        boolean[] visited = new boolean[customerCount + 1]; // 标记客户是否已访问
        Arrays.fill(visited, false);

        for (int i = 0; i < vehicleCount; i++) {
            double currentLoad = 0;
            int currentLocation = 0; // 仓库位置

            while (true) {
                double nearestDistance = Double.MAX_VALUE;
                int nearestCustomer = -1;
                for (int j = 1; j < demandArray.length; j++) {
                    if (!visited[j] && currentLoad + demandArray[j] <= capacity && distanceMatrix[currentLocation][j] < nearestDistance) {
                        nearestDistance = distanceMatrix[currentLocation][j];
                        nearestCustomer = j;
                    }
                }

                if (nearestCustomer == -1) break; // 没有可访问的客户

                // 更新解
                solution[index++] = nearestCustomer; // 添加客户
                visited[nearestCustomer] = true; // 标记客户为已访问
                totalDistance += nearestDistance; // 增加到最近客户的距离
                currentLoad += demandArray[nearestCustomer]; // 更新当前负载
                currentLocation = nearestCustomer; // 更新当前位置
            }

            totalDistance += distanceMatrix[currentLocation][0]; // 返回仓库的距离
            solution[index++] = 0; // 返回仓库
        }

        // 将总路径长度放到最后
        solution[index] = totalDistance; // 设置总路径长度
        return solution;
    }
    //节约算法（Savings Algorithm）来生成一个初始解
    public static double[] savingsAlgorithm(double[][] distanceMatrix, double[] demandArray) {
        ArrayList<int[]> savingsList = new ArrayList<>();

        // 计算节约值
        for (int i = 1; i <= customerCount; i++) {
            for (int j = i + 1; j <= customerCount; j++) {
                double savings = distanceMatrix[0][i] + distanceMatrix[0][j] - distanceMatrix[i][j];
                savingsList.add(new int[]{i, j, (int) savings});
            }
        }

        // 按节约值排序
        savingsList.sort((a, b) -> b[2] - a[2]); // 按降序排序

        boolean[] visited = new boolean[customerCount + 1];
        double[] solution = new double[vehicleCount + customerCount + 2];
        solution[0] = 0;
        int index = 1;
        double totalDistance = 0;

        for (int i = 0; i < vehicleCount; i++) {
            double currentLoad = 0; // 初始化当前载货量
            boolean vehicleHasVisited = false; // 标记当前车辆是否访问过客户
            int currentLocation = 0; // 当前车辆位置，初始化为仓库

            for (int[] saving : savingsList) {
                int customer1 = saving[0];
                int customer2 = saving[1];

                if (!visited[customer1] && !visited[customer2]) {
                    if (currentLoad + demandArray[customer1] + demandArray[customer2] <= capacity) {
                        // 将两个客户分配给当前车辆
                        visited[customer1] = true;
                        visited[customer2] = true;

                        solution[index++] = customer1;
                        solution[index++] = customer2;
                        totalDistance += distanceMatrix[currentLocation][customer1] +
                                distanceMatrix[customer1][customer2] +
                                distanceMatrix[currentLocation][customer2];
                        currentLoad += demandArray[customer1] + demandArray[customer2];
                        vehicleHasVisited = true;
                        currentLocation = customer2; // 更新当前位置为customer2
                    }
                } else if (!visited[customer1] && visited[customer2]) {
                    if (currentLoad + demandArray[customer1] <= capacity) {
                        // 将customer1分配给customer2的路径
                        visited[customer1] = true;
                        solution[index++] = customer1;
                        totalDistance += distanceMatrix[currentLocation][customer1];
                        currentLoad += demandArray[customer1];
                        vehicleHasVisited = true;
                        currentLocation = customer1; // 更新当前位置为customer1
                    }
                } else if (visited[customer1] && !visited[customer2]) {
                    if (currentLoad + demandArray[customer2] <= capacity) {
                        // 将customer2分配给customer1的路径
                        visited[customer2] = true;
                        solution[index++] = customer2;
                        totalDistance += distanceMatrix[currentLocation][customer2];
                        currentLoad += demandArray[customer2];
                        vehicleHasVisited = true;
                        currentLocation = customer2; // 更新当前位置为customer2
                    }
                }
            }

            // 如果车辆有访问过客户，返回仓库
            if (vehicleHasVisited) {
                solution[index++] = 0; // 返回仓库
                totalDistance += distanceMatrix[currentLocation][0]; // 加上返回仓库的距离
            }
        }
        solution[solution.length-1] = totalDistance; // 设置总路径长度
        return solution;
    }
    //最近插入法
    public static double[] nearestInsertionAlgorithm(double[][] distanceMatrix, double[] demandArray) {
        boolean[] visited = new boolean[customerCount + 1]; // 记录每个客户是否被访问过
        double[] solution = new double[vehicleCount + customerCount + 2];
        solution[0] = 0; // 仓库为起点
        int index = 1; // solution的索引
        double totalDistance = 0; // 总行驶距离

        for (int i = 0; i < vehicleCount; i++) {
            double currentLoad = 0; // 当前车辆的载重
            int currentLocation = 0; // 当前车辆的初始位置为仓库
            boolean vehicleHasVisited = false; // 标记车辆是否访问过客户
            int startIndex = index; // 记录当前车辆的起始索引

            // 初始化路径，选择最近的客户作为起点
            int firstCustomer = -1;
            double nearestDistance = Double.MAX_VALUE;
            for (int customer = 1; customer <= customerCount; customer++) {
                if (!visited[customer] && demandArray[customer] <= capacity) {
                    double distance = distanceMatrix[0][customer];
                    if (distance < nearestDistance) {
                        nearestDistance = distance;
                        firstCustomer = customer;
                    }
                }
            }

            // 没有可行的客户时结束
            if (firstCustomer == -1) {
                break;
            }

            // 将第一个客户加入路径
            visited[firstCustomer] = true;
            solution[index++] = firstCustomer;
            totalDistance += nearestDistance;
            currentLoad += demandArray[firstCustomer];
            currentLocation = firstCustomer;
            vehicleHasVisited = true;

            // 重复插入客户，直到无法插入更多客户
            while (true) {
                int bestCustomer = -1;
                double bestInsertionCost = Double.MAX_VALUE;
                int bestInsertPosition = -1;

                // 寻找最优插入的客户和位置
                for (int customer = 1; customer <= customerCount; customer++) {
                    if (!visited[customer] && currentLoad + demandArray[customer] <= capacity) {
                        // 仅在当前车辆的有效插入区间内查找
                        for (int j = startIndex; j < index; j++) {
                            int routePoint1 = (int) solution[j]; // 路径中的一个点
                            int routePoint2 = (j + 1 < index) ? (int) solution[j + 1] : 0; // 下一个点或仓库

                            // 计算插入到routePoint1和routePoint2之间的成本
                            double insertionCost = distanceMatrix[routePoint1][customer] + distanceMatrix[customer][routePoint2] - distanceMatrix[routePoint1][routePoint2];

                            // 寻找最佳插入点
                            if (insertionCost < bestInsertionCost) {
                                bestInsertionCost = insertionCost;
                                bestCustomer = customer;
                                bestInsertPosition = j + 1; // 插入位置
                            }
                        }
                    }
                }

                // 如果没有可以插入的客户，结束此路径
                if (bestCustomer == -1) {
                    break;
                }

                // 插入客户到最佳位置
                visited[bestCustomer] = true;
                currentLoad += demandArray[bestCustomer];
                totalDistance += bestInsertionCost;

                // 将客户插入到solution数组中
                if (index < solution.length - 1) {
                    System.arraycopy(solution, bestInsertPosition, solution, bestInsertPosition + 1, index - bestInsertPosition);
                    solution[bestInsertPosition] = bestCustomer;
                    index++;
                }
            }

            // 返回仓库
            if (vehicleHasVisited) {
                solution[index++] = 0; // 返回仓库
                totalDistance += distanceMatrix[currentLocation][0];
            }
        }

        solution[solution.length - 1] = totalDistance; // 设置总路径长度在数组最后一位
        return solution;
    }

    //检查解的容量
    public static boolean checkCapacity(double[] solution, double[] demandArray) {
        double currentLoad = 0;

        // 忽略solution数组的最后一个元素，它是总路径长度
        for (int i = 1; i < solution.length - 1; i++) {  // 从1开始，跳过最开始的仓库0
            int customer = (int) solution[i];

            // 如果遇到仓库（编号为0），表示一辆车的路径结束，检查载重
            if (customer == 0) {
                if (currentLoad > capacity) {
                    return false; // 超过容量限制，返回false
                }
                currentLoad = 0; // 重置当前载重，准备下一辆车
            } else {
                currentLoad += demandArray[customer]; // 累加客户需求量
            }
        }

        // 检查最后一辆车的载重
        return currentLoad <= capacity;
    }

    public static boolean checkAllSolutionsCapacity(List<double[]> solutions, double[] demandArray) {
        for (double[] solution : solutions) {
            // 对每个解进行容量检查
            if (!checkCapacity(solution, demandArray)) {
                return false; // 如果有任何一个解不满足容量限制，返回false
            }
        }
        return true; // 所有解都满足容量限制，返回true
    }

    //随机生成一个解
    public static double[] generateRandomSolution(double[][] distanceMatrix, double[] demandArray) {
        boolean[] visited = new boolean[customerCount + 1]; // 记录每个客户是否被访问过
        double[] solution = new double[vehicleCount + customerCount + 2]; // 解数组
        solution[0] = 0; // 仓库为起点
        int index = 1; // solution的索引
        Random random = new Random();
        double totalDistance = 0; // 用于计算总路径长度

        // 创建客户索引并打乱顺序
        Integer[] customers = new Integer[customerCount];
        for (int i = 0; i < customerCount; i++) {
            customers[i] = i + 1; // 客户编号从1开始
        }
        Collections.shuffle(Arrays.asList(customers), random); // 打乱客户顺序

        for (int i = 0; i < vehicleCount; i++) {
            double currentLoad = 0; // 当前车辆的载重
            int currentLocation = 0; // 当前车辆的初始位置为仓库
            boolean vehicleHasVisited = false; // 标记车辆是否访问过客户

            for (Integer customer : customers) {
                // 检查是否可以访问该客户
                if (!visited[customer] && currentLoad + demandArray[customer] <= capacity) {
                    // 计算从当前位置到该客户的距离并累加
                    totalDistance += distanceMatrix[currentLocation][customer]; // 累加距离
                    visited[customer] = true; // 标记客户为已访问
                    solution[index++] = customer; // 将客户加入路径
                    currentLoad += demandArray[customer]; // 更新当前负载
                    currentLocation = customer; // 更新当前位置
                    vehicleHasVisited = true;
                }

                // 如果车辆的负载已满，结束当前路径
                if (currentLoad >= capacity) {
                    break;
                }
            }

            // 返回仓库，并计算返回的距离
            if (vehicleHasVisited) {
                totalDistance += distanceMatrix[currentLocation][0]; // 返回仓库的距离
                solution[index++] = 0; // 返回仓库
            }
        }

        // 将总路径长度放到最后
        solution[solution.length - 1] = totalDistance; // 设置总路径长度
        return solution;
    }

    //初始化解集
    public static ArrayList<double[]> initialSolutions(double[][] distanceMatrix, double[] demandArray){
        ArrayList<double[]> solutions=new ArrayList<>();
        solutions.add(greedyAlgorithm(distanceMatrix,demandArray));
        solutions.add(nearestInsertionAlgorithm(distanceMatrix,demandArray));
        solutions.add(savingsAlgorithm(distanceMatrix,demandArray));
        for(int i=0;i<initialSolutionCount-3;i++){
            solutions.add(generateRandomSolution(distanceMatrix,demandArray));
        }
        return solutions;
    }

    //对解集按照总路径长度排序
    public static  void sortSolutionsByTotalLength(List<double[]> solutions){
        Collections.sort(solutions, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return Double.compare(o1[o1.length-1],o2[o2.length-1]);
            }
        });
    }

    //遗传算法(且确保满足容量限制)
    public static void genetic(ArrayList<double[]> solutions,double[]demandArray,double[][]distanceMatrix,double[][]points){
        int halfSize = initialSolutionCount / 2;
        int differenceCount=(int)(customerCount*0.3);//允许移动的范围
        int maxIndex=solutions.get(0).length-3;
        // 只保留前一半的解
        solutions.subList(halfSize, solutions.size()).clear();
        while(solutions.size()<initialSolutionCount){
            int length=solutions.get(0).length;
            double[] offspring1,offspring2;
            //顺序交叉
            //double p1=0;
            double p1=1.0-fe/FE;
            double r1=r.nextDouble();
            if(r1<=p1){
                double[][] offspring=getCrossoverOffspring(solutions);
                offspring1=offspring[0];
                offspring2=offspring[1];
            }else{
                int rand1=r.nextInt(initialSolutionCount/2);
                int rand2;
                do{
                    rand2=r.nextInt(initialSolutionCount/2);
                }while(rand2==rand1);
                offspring1=Arrays.copyOf(solutions.get(rand1),length);
                offspring2=Arrays.copyOf(solutions.get(rand2),length);
            }

            fe+=2;
/*            //交换片段变异
            offspring1=segmentSwapMutate(offspring1);
            offspring2=segmentSwapMutate(offspring2);*/
/*            //翻转一个片段
            offspring1=flipSegment(offspring1);
            offspring2=flipSegment(offspring2);*/
                //交换一位
            offspring1=swapMutate(offspring1,1);
            offspring2=swapMutate(offspring2,1);


/*            //挪动元素的位置
            offspring1=moveElement(offspring1,1);
            offspring2=moveElement(offspring2,1);*/
            /*//移动0的位置
            offspring1=moveZeroMutate(offspring1);
            offspring2=moveZeroMutate(offspring2);*/
/*            int maxTryCount=3;
            int count=0;
            while(count<maxTryCount){
                if(!checkCapacity(offspring1,demandArray)){
                    offspring1=moveZeroMutate(offspring1);
                    count++;
                    fe++;
                }else{
                    solutions.add(offspring1);
                    break;
                }

            }
            count=0;
            while(count<maxTryCount){
                if(!checkCapacity(offspring2,demandArray)){
                    offspring2=moveZeroMutate(offspring2);
                    count++;
                    fe++;
                }else{
                    solutions.add(offspring2);
                    break;
                }

            }*/
            offspring1=clustering(offspring1,distanceMatrix,points);
            offspring2=clustering(offspring2,distanceMatrix,points);
            if(!checkCapacity(offspring1,demandArray)){
                offspring1=adjustAllocation(offspring1,demandArray,vehicleCount,distanceMatrix);
                fe++;
            }else {
                solutions.add(offspring1);
            }
            if(!checkCapacity(offspring2,demandArray)){
                offspring2=adjustAllocation(offspring2,demandArray,vehicleCount,distanceMatrix);
                fe++;
            }else {
                solutions.add(offspring2);
            }

        }

    }

    //对随机两辆车进行聚类，重新分配
    public static double[] clustering(double[]solution,double[][]distanceMatrix,double[][]points){
        List<List<Integer>> vehicleList=splitVehicle(solution);

        int time1=1;//选择两辆车的次数
        int cur1=0;
        int time2=10;//选中的两辆车聚类的次数
        int cur2=0;
        while(cur1<time1){
            /*int[] v=getV1AndV2(vehicleList,points);
            int v1=v[0],v2=v[1];*/
            int v1=r.nextInt(vehicleList.size()),v2;
            do{
                v2=r.nextInt(vehicleList.size());
            } while (v1==v2);

            List<Integer>vehicle1=vehicleList.get(v1);
            List<Integer>vehicle2=vehicleList.get(v2);
            while(cur2<time2){
                double[] center1=findClusterCenter(vehicle1,points);
                double[] center2=findClusterCenter(vehicle2,points);
                List<Integer> mergedList = new ArrayList<>(vehicle1);
                mergedList.addAll(vehicle2);
                vehicle1.clear();
                vehicle2.clear();
                for(int cust:mergedList){
                    double distance1=Math.pow(points[cust][0]-center1[0],2)+Math.pow(points[cust][1]-center1[1],2);
                    double distance2=Math.pow(points[cust][0]-center2[0],2)+Math.pow(points[cust][1]-center2[1],2);
                    if(distance1<=distance2){
                        vehicle1.add(cust);
                    }else vehicle2.add(cust);
                }
                cur2++;
            }
            cur1++;
        }


        while (vehicleList.size()<vehicleCount) vehicleList.add(new ArrayList<>());
        //先头尾都加上0，以方便后续计算
        for (int j = 0; j < vehicleList.size(); j++) {
            vehicleList.get(j).add(0,0);
            vehicleList.get(j).add(vehicleList.get(j).size(),0);
        }
        double[]newSolution=transferToSolution(vehicleList,solution,distanceMatrix);
        return newSolution;


    }

    //找最近的两辆车
    public static int[] getV1AndV2(List<List<Integer>> vehicleList,double[][]points){
        int n=vehicleList.size();
        double[][] distance=new double[n][n];
        double[][] centerList=new double[n][2];
        for (int i = 0; i < n; i++) {
            List<Integer> vehicle=vehicleList.get(i);
            double[] centerI=new double[2];
            for (int j = 0; j < vehicle.size(); j++) {
                centerI[0]+=points[vehicle.get(j)][0];
                centerI[1]+=points[vehicle.get(j)][1];
            }
            centerI[0]/=(double) vehicle.size();
            centerI[1]/=(double) vehicle.size();
            centerList[i]=centerI;
        }
      /*  double minDistanceSquared = Double.MAX_VALUE;
        int index1 = -1, index2 = -1;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double distanceSquared = Math.pow(centerList[i][0] - centerList[j][0], 2) +
                        Math.pow(centerList[i][1] - centerList[j][1], 2);
                if (distanceSquared < minDistanceSquared) {
                    minDistanceSquared = distanceSquared;
                    index1 = i;
                    index2 = j;
                }
            }
        }*/
        // 找到距离最远的两个中心点的索引
        double maxDistanceSquared = 0;
        int index1 = -1, index2 = -1;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double distanceSquared = Math.pow(centerList[i][0] - centerList[j][0], 2) +
                        Math.pow(centerList[i][1] - centerList[j][1], 2);
                if (distanceSquared > maxDistanceSquared) {
                    maxDistanceSquared = distanceSquared;
                    index1 = i;
                    index2 = j;
                }
            }
        }



        // 返回最近两个中心点的索引
        return new int[]{index1, index2};


    }
    //找聚类中心
    public static double[]findClusterCenter(List<Integer>vehicle,double[][]points){
        double[] center=new double[2];
        for(int cust:vehicle){
            center[0]+=points[cust][0];
            center[1]+=points[cust][1];
        }
        double size=vehicle.size();
        if(size!=0){
            center[0]/=size;
            center[1]/=size;
        }

        return center;
    }

    //不满足容量限制的解，将每辆车看做一类，对于超出容量限制的车辆，拿出离其他客户距离最远的点来放入离他最近的车辆，如果不满足容量限制就次近
    //adjustAllocation直接计算好超出容量限制的车辆里的客户与其他客户的距离，不随着移动来重新距离，节约时间
    public static double[]  adjustAllocation(double[]solution,double[]demandArray,int vehicleCount,double[][]distanceMatrix){
        List<List<Integer>> vehicleList=splitVehicle(solution);
        while (vehicleList.size()<vehicleCount) vehicleList.add(new ArrayList<>());
        List<Double> vehicleDemandList=new ArrayList<>();
        List<Integer> exceededVehicle=new ArrayList<>();//超出容量的车辆
        //每辆车的需求量
        for(int i=0;i<vehicleList.size();i++){
            List<Integer> vehicle=vehicleList.get(i);
            double sum=0;
            for(int cust:vehicle){
               sum+=demandArray[cust];
            }
            if (sum>capacity) exceededVehicle.add(i);
            vehicleDemandList.add(sum);
        }

        //先选出要拿出来的客户
        List<Integer>removeCustList=new ArrayList<>();
        for(int i=0;i<exceededVehicle.size();i++){//对每一个超出容量限制的车辆
            int vehicleInx=exceededVehicle.get(i);
            List<Integer> vehicle=vehicleList.get(vehicleInx);
            //计算超出容量限制的车辆里的每个客户到其他客户的距离
            List<Double> vehicleDistance=new ArrayList<>();
            for (int j = 0; j < vehicle.size(); j++) {
                double distanceSum=0;
                int cust=vehicle.get(j);
                for (int k = 0; k < vehicle.size(); k++) {
                    distanceSum+=distanceMatrix[cust][vehicle.get(k)];
                }
                vehicleDistance.add(distanceSum);
            }
            while(vehicleDemandList.get(vehicleInx)>capacity){
                double maxDistanceSum=0;
                int selectCust=-1;
                int selectCustInx=-1;
                for (int j = 0; j < vehicleDistance.size(); j++) {
                    if(vehicleDistance.get(j)>maxDistanceSum){
                        maxDistanceSum=vehicleDistance.get(j);
                        selectCust=vehicle.get(j);
                        selectCustInx=j;
                    }
                }
                vehicle.remove(selectCustInx);
                removeCustList.add(selectCust);
                vehicleDemandList.set(vehicleInx,vehicleDemandList.get(vehicleInx)-demandArray[selectCust]);
                vehicleDistance.remove(selectCustInx);
            }

        }
        //对需要重新插入的客户进行从大到小排序，以防止需求较大的客户插不进去的情况
        Collections.sort(removeCustList,(a,b)->Double.compare(demandArray[b],demandArray[a]));

        //现在将客户插入其他车辆中
        for (int i = 0; i < removeCustList.size(); i++) {
            int cust=removeCustList.get(i);
            double minAverge=Double.MAX_VALUE;
            int vehicleInx=-1;
            for (int j = 0; j < vehicleDemandList.size(); j++) {
                if(demandArray[cust]+vehicleDemandList.get(j)>capacity) continue;
                double avergeDistance=0;
                List<Integer> vehicle=vehicleList.get(j);
                for(int customer:vehicle){
                    avergeDistance+=distanceMatrix[cust][customer];
                }
                if(vehicle.size()!=0)
                    avergeDistance/=(double) vehicle.size();
                if(avergeDistance<minAverge){
                    minAverge=avergeDistance;
                    vehicleInx=j;
                }
            }
            //当前cust装不下去其他车辆，就随机选一辆车，一直往外拿客户，直到能装下当前客户
            if(vehicleInx==-1){
                int randomVehicleInx=r.nextInt(vehicleCount);
                List<Integer> randomVehicle=vehicleList.get(randomVehicleInx);
                while(vehicleDemandList.get(randomVehicleInx)>capacity-demandArray[cust]){
                    int removeCust=randomVehicle.remove(0);
                    vehicleDemandList.set(randomVehicleInx,vehicleDemandList.get(randomVehicleInx)-demandArray[removeCust]);
                    removeCustList.add(removeCust);
                }
                vehicleInx=randomVehicleInx;
            }
            List<Integer>vehicle=vehicleList.get(vehicleInx);
            vehicle.add(cust);
            vehicleDemandList.set(vehicleInx,demandArray[cust]+vehicleDemandList.get(vehicleInx));
        }
        //先头尾都加上0，以方便后续计算
        for (int j = 0; j < vehicleList.size(); j++) {
            vehicleList.get(j).add(0,0);
            vehicleList.get(j).add(vehicleList.get(j).size(),0);
        }
        double[]newSolution=transferToSolution(vehicleList,solution,distanceMatrix);
        return newSolution;



    }

    //adjustAllocation1每拿出一个客户，都重新计算距离
    public static double[]  adjustAllocation1(double[]solution,double[]demandArray,int vehicleCount,double[][]distanceMatrix){
        List<List<Integer>> vehicleList=splitVehicle(solution);
        while (vehicleList.size()<vehicleCount) vehicleList.add(new ArrayList<>());
        List<Double> vehicleDemandList=new ArrayList<>();
        List<Integer> exceededVehicle=new ArrayList<>();//超出容量的车辆
        //每辆车的需求量
        for(int i=0;i<vehicleList.size();i++){
            List<Integer> vehicle=vehicleList.get(i);
            double sum=0;
            for(int cust:vehicle){
                sum+=demandArray[cust];
            }
            if (sum>capacity) exceededVehicle.add(i);
            vehicleDemandList.add(sum);
        }

        //先选出要拿出来的客户
        List<Integer>removeCustList=new ArrayList<>();
        for(int i=0;i<exceededVehicle.size();i++){//对每一个超出容量限制的车辆
            int vehicleInx=exceededVehicle.get(i);
            while (vehicleDemandList.get(vehicleInx)>capacity){//一直往外客户，直到满足限制
                double maxDistanceSum=0;
                int selectCust=-1;
                int selectCustInx=-1;
                //选出离得最远的点，直到满足容量限制
                List<Integer> vehicle=vehicleList.get(vehicleInx);
                for (int j = 0; j < vehicle.size(); j++) {
                    double distanceSum=0;
                    int cust=vehicle.get(j);
                    for (int k = 0; k < vehicle.size(); k++) {
                        if(j==k)continue;
                        distanceSum+=distanceMatrix[cust][vehicle.get(k)];
                    }
                    if(distanceSum>maxDistanceSum){
                        maxDistanceSum=distanceSum;
                        selectCust=cust;
                        selectCustInx=j;
                    }
                }
                removeCustList.add(selectCust);
                vehicle.remove(selectCustInx);
                vehicleDemandList.set(vehicleInx,vehicleDemandList.get(vehicleInx)-demandArray[selectCust]);
            }
        }
        //现在将客户插入其他车辆中
        for (int i = 0; i < removeCustList.size(); i++) {
            int cust=removeCustList.get(i);
            double minAverge=Double.MAX_VALUE;
            int vehicleInx=-1;
            for (int j = 0; j < vehicleDemandList.size(); j++) {
                if(demandArray[cust]+vehicleDemandList.get(j)>capacity) continue;
                double avergeDistance=0;
                List<Integer> vehicle=vehicleList.get(j);
                for(int customer:vehicle){
                    avergeDistance+=distanceMatrix[cust][customer];
                }
                if(vehicle.size()!=0)
                    avergeDistance/=(double) vehicle.size();
                if(avergeDistance<minAverge){
                    minAverge=avergeDistance;
                    vehicleInx=j;
                }
            }
            List<Integer>vehicle=vehicleList.get(vehicleInx);
            vehicle.add(cust);
            vehicleDemandList.set(vehicleInx,demandArray[cust]+vehicleDemandList.get(vehicleInx));


        }
        //先头尾都加上0，以方便后续计算
        for (int j = 0; j < vehicleList.size(); j++) {
            vehicleList.get(j).add(0,0);
            vehicleList.get(j).add(vehicleList.get(j).size(),0);
        }
        double[]newSolution=transferToSolution(vehicleList,solution,distanceMatrix);
        return newSolution;


    }
    //交叉得到子代
    public static double[][] getCrossoverOffspring(List<double[]>solutions) {
        Random r=new Random();
        // 从 solutions 中随机选择两个父母
        int parentIndex1 = r.nextInt(solutions.size());
        int parentIndex2;
        do {
            parentIndex2 = r.nextInt(solutions.size());
        } while (parentIndex1 == parentIndex2);  // 确保两个父母不同

        double[] parent1 = solutions.get(parentIndex1);
        double[] parent2 = solutions.get(parentIndex2);
        int length = parent1.length;

        // 确定交叉片段范围（跳过第一个和最后两个元素）
        int start = 1 + r.nextInt(length - 3); // start 范围为 [1, length - 3]
        int end = start + r.nextInt(length - start - 2) + 1; // end 范围为 [start, length - 2]

        double[]offspring1=crossover(start,end,parent1,parent2);
        double[]offspring2=crossover(start,end,parent2,parent1);

        return new double[][]{offspring1, offspring2}; // 返回两个子代
    }
    //交叉操作
    public static double[] crossover(int start,int end,double[]parent1,double[]parent2){
        double[] offspring=new double[parent1.length];
        Set<Double> set=new HashSet<>();
        // 1. 生成子代1
        for(int i=start;i<=end;i++){
            if(parent1[i]!=0){
                offspring[i]=parent1[i];
                set.add(parent1[i]);
            }
        }

        // 填充剩余的空位
        int idx = 1; // 用于遍历 parent2
        for(int i=1;i<offspring.length-2;i++){
            if(offspring[i]==0){
                while(set.contains(parent2[idx])) idx++;
                offspring[i]=parent2[idx++];
            }
        }
        return offspring;
    }

    public static double[] moveZeroMutate(double[] originSolution){
        int maxIndex=originSolution.length-3;
        int differenceCount=(int)(customerCount*0.4);//允许移动的范围，最多前后移动customerCount*0.3个客户位置
        double[] mutatedSolution = originSolution.clone();
        // 随机选择一个索引范围，确保不包括第一个和最后两个元素
        int a = r.nextInt(maxIndex) + 1; // start 范围为 [1, selectedSolution.length - 3]
        int b = r.nextInt(maxIndex) + 1;
        int start=Math.min(a,b);
        int end=Math.max(a,b);
        //为了照顾一下最后一辆车
        double p=r.nextDouble();
        if(p<0.7 && mutatedSolution[maxIndex]==0) end=maxIndex;

        ArrayList<Integer> zeroIndices = new ArrayList<>();
        for (int j = start; j <= end; j++) {
            if (mutatedSolution[j] == 0.0) {
                zeroIndices.add(j);
            }
        }
        // 如果找到0的位置，进行交换
        for(int j=0;j<zeroIndices.size();j++){
            int zeroIndex=zeroIndices.get(j);
            int randomIndexToSwap;
            do {
                randomIndexToSwap = r.nextInt(differenceCount*2) + zeroIndex-differenceCount;

            } while ( randomIndexToSwap>maxIndex || randomIndexToSwap<1); // 确保选择的索引合法

            double temp=mutatedSolution[zeroIndex];
            mutatedSolution[zeroIndex]=mutatedSolution[randomIndexToSwap];
            mutatedSolution[randomIndexToSwap]=temp;
        }
        return mutatedSolution;
    }

    //交换变异
    public static double[] swapMutate(double[] originSolution,int swapCount) {
        // 创建原解的副本，防止修改原解
        double[] mutatedSolution = originSolution.clone();
        int count=0;
        while(count<swapCount){
            // 随机选择两个不同的索引进行交换
            int index1 = r.nextInt(mutatedSolution.length-3)+1;
            int index2;
            do {
                index2 = r.nextInt(mutatedSolution.length-3)+1;
            } while (index2 == index1); // 确保 index2 不等于 index1

            // 交换这两个索引的值
            double temp = mutatedSolution[index1];
            mutatedSolution[index1] = mutatedSolution[index2];
            mutatedSolution[index2] = temp;
            count++;
        }


        return mutatedSolution; // 返回变异后的解
    }

    //随机选一个元素挪动位置
    public static double[] moveElement(double[] originSolution, int moveCount) {
        double[] mutatedSolution = originSolution.clone();
        int count = 0;

        while (count < moveCount) {
            // 随机选择要移动的元素的索引
            int index1 = r.nextInt(mutatedSolution.length - 3) + 1; // index1 范围从 1 到 length-3
            int index2;

            // 随机选择新位置的索引，确保不与 index1 相同
            do {
                index2 = r.nextInt(mutatedSolution.length - 3) + 1;
            } while (index2 == index1); // 确保 index2 不等于 index1

            // 移动元素
            double temp = mutatedSolution[index1]; // 保存要移动的元素
            // 向后移动 index1 到 index2 之间的元素
            if (index2 < index1) {
                System.arraycopy(mutatedSolution, index2, mutatedSolution, index2 + 1, index1 - index2);
            } else {
                System.arraycopy(mutatedSolution, index1 + 1, mutatedSolution, index1, index2 - index1);
            }
            mutatedSolution[index2] = temp; // 将元素放到新位置

            count++;
        }

        return mutatedSolution; // 返回变异后的解
    }

    //交换片段
    public static double[] segmentSwapMutate(double[] originSolution) {
        // 创建原解的副本，防止修改原解
        double[] mutatedSolution = originSolution.clone();
        int length = mutatedSolution.length;

        // 随机选择4个不同的索引，避开第一个和最后两个元素
        Set<Integer> indices = new HashSet<>();
        while (indices.size() < 4) {
            indices.add(r.nextInt(length - 3) + 1); // 从1到length-3
        }

        // 将索引转换为数组并排序
        Integer[] indexArray = indices.toArray(new Integer[0]);
        Arrays.sort(indexArray);

        // 定义片段1和片段2
        int start1 = indexArray[0]; // 最小的索引
        int end1 = indexArray[1];   // 第二小的索引
        int start2 = indexArray[2]; // 第三个索引
        int end2 = indexArray[3];   // 最大的索引

        // 创建一个新的数组来存储变异后的解
        double[] newSolution = new double[length];

        // 复制片段1前面的元素
        System.arraycopy(mutatedSolution, 0, newSolution, 0, start1);

        // 复制片段2的元素
        System.arraycopy(mutatedSolution, start2, newSolution, start1, end2 - start2 + 1);

        // 复制片段1之后且片段2之前的元素
        System.arraycopy(mutatedSolution, end1 + 1, newSolution, start1 + (end2 - start2 + 1), start2 - end1 - 1);

        // 复制片段1的元素
        System.arraycopy(mutatedSolution, start1, newSolution, start1 + (end2 - start2 + 1) + (start2 - end1 - 1), end1 - start1 + 1);

        // 复制片段2之后的元素
        System.arraycopy(mutatedSolution, end2 + 1, newSolution, start1 + (end2 - start2 + 1) + (start2 - end1 - 1) + (end1 - start1 + 1), length - end2 - 1);

        return newSolution; // 返回变异后的解
    }

    //选一个片段翻转
    public static double[] flipSegment(double[] originSolution) {
        // 创建原解的副本，防止修改原解
        double[] mutatedSolution = originSolution.clone();

        // 获取原解的长度
        int length = originSolution.length;
        int maxIndex = length - 3; // 片段选择的最大索引

        // 随机选择一个片段的起始位置和结束位置，确保不包括第一个和最后两个元素
        int start = r.nextInt(maxIndex) + 1; // start 范围为 [1, length - 3]
        int end = start + r.nextInt(maxIndex - start + 1); // end 范围为 [start, length - 3]

        // 进行片段翻转
        while (start < end) {
            double temp = mutatedSolution[start];
            mutatedSolution[start] = mutatedSolution[end];
            mutatedSolution[end] = temp;
            start++;
            end--;
        }

        return mutatedSolution; // 返回变异后的解
    }


    //初始化信息素矩阵
    public static double[][] initialPherome(double deltaTau){
        double initialValue=deltaTau;
        double[][] pheromeMatrix=new double[customerCount+1][customerCount+1];
        for(double[] row:pheromeMatrix){
            Arrays.fill(row,initialValue);
        }
        for(int i=0;i<customerCount+1;i++)
            pheromeMatrix[i][i]=0;
        return pheromeMatrix;
    }

    //蚁群排序
    public static void ACSSort(List<double[]> solutions,double[][]pheromoneMatrix,double[][] distanceMatrix,List<double[]>archive){
        for(int i=0;i<solutions.size();i++){
            double[] solution=solutions.get(i);
            List<List<Integer>> vehicleList=splitVehicle(solution);
            boolean[] visited=new boolean[customerCount+1];//第一个放仓库
            visited[0]=true;
            double totalLength=0;
            int index=1;
            for(int j=0;j<vehicleList.size();j++){
                List<Integer> vehicle=vehicleList.get(j);
                int currentCustomer=0;
                while(!isAllCustomerOfVehicleVisited(visited,vehicle)){
                   int nextCustomer=getNextCustomer(visited,vehicle,pheromoneMatrix,distanceMatrix,currentCustomer);
                   solution[index++]=nextCustomer;
                   totalLength+=distanceMatrix[currentCustomer][nextCustomer];
                   visited[nextCustomer]=true;
                   //局部信息变化
                    pheromoneMatrix[currentCustomer][nextCustomer] = (1 - rho) * pheromoneMatrix[currentCustomer][nextCustomer] + rho * deltaTau;
                   currentCustomer=nextCustomer;
                }
                //回到仓库
                solution[index++]=0;
                totalLength+=distanceMatrix[currentCustomer][0];
            }
            //以防没有用完所有车辆
            while(index<solution.length-1){
                solution[index++]=0;
            }
            solution[solution.length-1]=totalLength;
        }
        //全局信息更新
        sortSolutionsByTotalLength(solutions);
        double[]solution0=archive.get(0);
        double value0=solution0[solution0.length-1];
        if(value0<bestValue){
            bestSolution=Arrays.copyOf(solution0,solution0.length);
            bestValue=value0;
        }
        for(int i=0;i<bestSolution.length-2;i++){
            int front=(int)bestSolution[i];
            int behind=(int)bestSolution[i+1];
            pheromoneMatrix[front][behind]+=alpha*(1.0/bestValue);
        }




    }

    //局部搜索1,拿掉每一辆车的一个客户，使得减少的路程长度最多,将客户重新插到车辆中,使得增加的路程最少
    public static void localSearch(List<double[]> solutions,double[][]distanceMatrix,double[]demandArray){
        for (int i = 0; i < solutions.size(); i++) {
            double[] solution=solutions.get(i);
            List<List<Integer>> vehicleList=splitVehicle(solution);
            //先头尾都加上0，以方便后续计算
            for (int j = 0; j < vehicleList.size(); j++) {
                vehicleList.get(j).add(0,0);
                vehicleList.get(j).add(vehicleList.get(j).size(),0);
            }
            //先得到要移除的客户列表
            List<Integer> removeCustList=getRemoveCustList(vehicleList,distanceMatrix);
            //将客户插入到车辆中
            insertCust(vehicleList,removeCustList,distanceMatrix,demandArray);
            //将修改过的车辆的客户列表转变回解
            solution=transferToSolution(vehicleList,solution,distanceMatrix);

        }
    }
    //局部搜索2，将一辆车（路径最长的）里面的客户全部重新分配给其他车辆，然后再进行其他调整以满足容量限制
    public static ArrayList<double[]> localSearch2(List<double[]> solutions,double[][]distanceMatrix,double[]demandArray){
        ArrayList<double[]> newSolutions=new ArrayList<>();
        for (int i = 0; i < solutions.size(); i++) {
            double[] solution=Arrays.copyOf(solutions.get(i),solutions.get(0).length);
            List<List<Integer>> vehicleList=splitVehicle(solution);
            /*//随机选一辆车，将里面的所有客户重新分配

            int randV=r.nextInt(vehicleList.size());*/
            int randV=getTheLengthMaxVehicle(distanceMatrix,vehicleList);
            List<Integer> removeCustList=new ArrayList<>(vehicleList.get(randV));
            vehicleList.get(randV).clear();
            while(vehicleList.size()<vehicleCount) vehicleList.add(new ArrayList<>());
            //将客户插入到车辆中
            //先头尾都加上0，以方便后续计算
            for (int j = 0; j < vehicleList.size(); j++) {
                vehicleList.get(j).add(0,0);
                vehicleList.get(j).add(vehicleList.get(j).size(),0);
            }
            insertCust(vehicleList,removeCustList,distanceMatrix,demandArray);
            //将修改过的车辆的客户列表转变回解
            solution=transferToSolution(vehicleList,solution,distanceMatrix);
            /*if(!checkCapacity(solution,demandArray)){
                solution=adjustAllocation(solution,demandArray,vehicleCount,distanceMatrix);
                fe++;
            }*/

           newSolutions.add(Arrays.copyOf(solution,solution.length));
        }
        return newSolutions;
    }

    //选出路径最长的一辆车
    public static int getTheLengthMaxVehicle(double[][]distanceMatrix,List<List<Integer>> vehicleList){
        List<List<Integer>> newVehicleList = new ArrayList<>();
        for (List<Integer> vehicle : vehicleList) {
            newVehicleList.add(new ArrayList<>(vehicle)); // 深拷贝内层列表
        }
        double maxLength=0;
        int v=-1;
        for (int i = 0; i < newVehicleList.size(); i++) {
            //先给每辆车加前后仓库方便计算
            List<Integer> vehicle=newVehicleList.get(i);
            vehicle.add(0,0);
            vehicle.add(vehicle.size(),0);

            double length=0;
            for (int j = 0; j < vehicle.size()-1; j++) {
                length+=distanceMatrix[vehicle.get(j)][vehicle.get(j+1)];
            }
            if(length>maxLength){
                maxLength=length;
                v=i;
            }
        }
        return v;
    }
    //拿掉每一辆车的一个客户，使得减少的路程长度最多
    public static List<Integer> getRemoveCustList(List<List<Integer>> vehicleList,double[][]distanceMatrix) {
        List<Integer> removeCustList = new ArrayList<>();
        for (int j = 0; j < vehicleList.size(); j++) {
            List<Integer> vehicle = vehicleList.get(j);
            if (vehicle.size() <=3) continue;//如果一辆车只有一个客户，那不用拿掉一个客户
            else {
                double maxDifference = 0;
                int cust = -1;
                for (int k = 2; k < vehicle.size() - 1; k++) {//不考虑拿掉第一个
                    double difference = distanceMatrix[vehicle.get(k - 1)][vehicle.get(k)]
                            + distanceMatrix[vehicle.get(k)][vehicle.get(k + 1)]
                            - distanceMatrix[vehicle.get(k - 1)][vehicle.get(k + 1)];
                    if (difference > maxDifference) {
                        cust = vehicle.get(k);
                        maxDifference = difference;
                    }
                }
                removeCustList.add(cust);
                vehicleList.get(j).remove(Integer.valueOf(cust));
            }

        }
        return removeCustList;

    }

    //将客户重新插到车辆中,使得增加的路程最少，考虑容量限制
    public static void insertCust(List<List<Integer>> vehicleList, List<Integer> removeCustList, double[][] distanceMatrix, double[] demandArray) {
        for (int k = 0; k < removeCustList.size(); k++) {
            int cust = removeCustList.get(k);
            double minDifference = Double.MAX_VALUE;
            int veInx = -1;
            int inInx = -1;

            for (int i = 0; i < vehicleList.size(); i++) {
                List<Integer> vehicle = vehicleList.get(i);
                double demandSum = 0;

                // 计算当前车辆的需求总和
                for (int j = 0; j < vehicle.size(); j++) {
                    demandSum += demandArray[vehicle.get(j)];
                }

                // 检查是否能容纳新客户
                if (demandSum + demandArray[cust] > capacity) continue;

                // 查找插入位置
                for (int j = 0; j < vehicle.size() - 1; j++) {
                    double difference = distanceMatrix[vehicle.get(j)][cust]
                            + distanceMatrix[cust][vehicle.get(j + 1)]
                            - distanceMatrix[vehicle.get(j)][vehicle.get(j + 1)];

                    if (difference < minDifference) {
                        minDifference = difference;
                        veInx = i;
                        inInx = j + 1;
                    }
                }
            }

            // 确保找到合适的车辆和插入位置
            if (veInx != -1 && inInx != -1) {
                vehicleList.get(veInx).add(inInx, cust);
            }
        }
    }

    //将客户重新插到车辆中,使得增加的路程最少，不考虑容量限制
    public static void insertCustRegardlessCapa(List<List<Integer>> vehicleList, List<Integer> removeCustList, double[][] distanceMatrix, double[] demandArray) {
        for (int k = 0; k < removeCustList.size(); k++) {
            int cust = removeCustList.get(k);
            double minDifference = Double.MAX_VALUE;
            int veInx = -1;
            int inInx = -1;

            for (int i = 0; i < vehicleList.size(); i++) {
                List<Integer> vehicle = vehicleList.get(i);
                double demandSum = 0;

                // 计算当前车辆的需求总和
                for (int j = 0; j < vehicle.size(); j++) {
                    demandSum += demandArray[vehicle.get(j)];
                }

                // 查找插入位置
                for (int j = 0; j < vehicle.size() - 1; j++) {
                    double difference = distanceMatrix[vehicle.get(j)][cust]
                            + distanceMatrix[cust][vehicle.get(j + 1)]
                            - distanceMatrix[vehicle.get(j)][vehicle.get(j + 1)];

                    if (difference < minDifference) {
                        minDifference = difference;
                        veInx = i;
                        inInx = j + 1;
                    }
                }
            }

            // 确保找到合适的车辆和插入位置
            if (veInx != -1 && inInx != -1) {
                vehicleList.get(veInx).add(inInx, cust);
            }
        }
    }


    //将车辆列表变为解
    public static double[] transferToSolution(List<List<Integer>> vehicleList, double[] solution, double[][] distanceMatrix) {
        double[] newSolution = new double[solution.length];
        int index = 1;

        for (List<Integer> vehicle : vehicleList) {
            for (int j = 1; j < vehicle.size() - 1; j++) {
                newSolution[index++] = vehicle.get(j);
            }
            newSolution[index++] = 0; // 添加回到起点的0
        }

        // 计算总路径长度
        double totalLength = 0;
        for (int i = 0; i < index - 1; i++) {
            totalLength += distanceMatrix[(int)newSolution[i]][(int)newSolution[i + 1]];
        }
        newSolution[newSolution.length - 1] = totalLength; // 最后一项为总路径长度

        return newSolution;
    }
    public static double cucalateTotalLength(double[]solution,double[][]distanceMatrix){
        double totalLength = 0;
        for (int i = 0; i < solution.length - 2; i++) {
            totalLength += distanceMatrix[(int)solution[i]][(int)solution[i + 1]];
        }
        return totalLength;
    }

    //得到每一辆车的客户序列
    public static List<List<Integer>> splitVehicle(double[]solution){
        List<List<Integer>> vehicleList=new ArrayList<>();
        List<Integer> vehicle=new ArrayList<>();
        for(int i=1;i<solution.length-1;i++){
            if(solution[i]==0){
                vehicleList.add(vehicle);
                vehicle=new ArrayList<>();
            }else{
                vehicle.add((int)solution[i]);
            }
        }
        return vehicleList;
    }

    //检查一辆车的所有客户是不是已经访问完
    public static boolean isAllCustomerOfVehicleVisited(boolean[]visited,List<Integer> vehicle){
        for(int customer:vehicle){
            if(!visited[customer]) return false;
        }
        return true;
    }

    //蚁群算法选择下一个访问的城市
    public static int getNextCustomer(boolean[] visited,List<Integer> vehicle, double[][]pheromoneMatrix, double[][] distanceMatrix, int currentCustomer) {
        Random r = new Random();
        int nextCust = -1;//初始化下一个访问城市
        double q = r.nextDouble();
        while (nextCust == -1) {
            if (q <= q0)
            {
                double max = 0;
                for (int k = 0; k < vehicle.size(); k++) {
                    int custk=vehicle.get(k);
                    if (!visited[custk]) {
                        double temp = 1 * 1e60 * pheromoneMatrix[currentCustomer][custk] * Math.pow((1.0 / distanceMatrix[currentCustomer][custk]), beta);
                        if (temp > max) {
                            max = temp;
                            nextCust = custk;
                        }
                    }
                }

            } else {
                double[] selectP = new double[vehicle.size()];//用于存放未访问城市的选择概率
                double fenmu = 0;
                for (int k = 0; k < vehicle.size(); k++) {
                    int custk = vehicle.get(k);
                    if (!visited[custk]) {
                        double fenzi = 1e30 * pheromoneMatrix[currentCustomer][custk] * Math.pow(1.0 / distanceMatrix[currentCustomer][custk], beta);
                        fenmu += fenzi;
                        selectP[k] = fenzi;
                    }
                }
                for (int k = 0; k < selectP.length; k++) {
                    selectP[k] /= fenmu;
                }
                double p = r.nextDouble();
                double sum = 0;
                for (int k = 0; k < selectP.length; k++) {
                    sum += selectP[k];
                    if (p <= sum) {
                        nextCust = vehicle.get(k);
                        break;
                    }
                }
            }
           /* if(nextCust==-1)
                System.out.println("-----------------------");
*/
        }
        return nextCust;
    }


    //更新Archive
    public static List<double[]> updateArchive(List<double[]> archive, List<double[]> solutions) {
        List<double[]> newArchive = new ArrayList<>(archive);
        newArchive.addAll(solutions);

        // 使用 Set 去重
        Set<String> uniqueSet = new HashSet<>();
        List<double[]> distinctArchive = new ArrayList<>();

        for (double[] solution : newArchive) {
            String key = Arrays.toString(solution); // 将数组转换为字符串
            if (uniqueSet.add(key)) { // 如果集合中没有该键，则添加
                distinctArchive.add(Arrays.copyOf(solution,solution.length));
            }
        }

        sortSolutionsByTotalLength(distinctArchive);
        double maxValue = bestValue * 1.02;
        int count = 0;
        int lengthIndex = distinctArchive.get(0).length - 1;

        while (count < distinctArchive.size() && distinctArchive.get(count)[lengthIndex] <=maxValue) {
            count++;
        }

        distinctArchive.subList(count, distinctArchive.size()).clear();
        return distinctArchive;
    }

    //打印archive
    public static void printArchive(List<double[]> archive) {
        for (double[] array : archive) {
            System.out.println(Arrays.toString(array)); // 打印每个 double[] 数组
        }
    }

    //保存解
    public static void saveArchive(List<double[]> archive, String saveFilePath) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(saveFilePath, true))) {
            for (int i = 0; i < 1; i++) {//我现在只写每一次运行的最优解
                double[] solution = archive.get(i);
                bw.write("length:" + solution[solution.length - 1]); // 修正为转换长度为整数
                bw.newLine();
                bw.write(0+",");
                for (int j = 1; j < solution.length - 2; j++) {
                    if((int)solution[j]==0){
                        bw.write(0+"\r\n");
                        bw.write(0+",");
                    }
                    else bw.write((int)solution[j]+",");
                }
                bw.write(0+"\r\n");
            }
            bw.newLine();
            bw.newLine();
        } catch (IOException e) {
            System.out.println("保存数据异常: " + e.getMessage()); // 输出异常信息以便调试
        }
    }

    //对Archive里面的第一个解进行局部搜索，使得更多类似的解
    public static List<double[]> archive0LocalSearch(List<double[]>archive,double[][]points,double[][]distanceMatrix){
        double[]solution=Arrays.copyOf(archive.get(0),archive.get(0).length);
        List<double[]> newSolutions=new ArrayList<>(archive);
        int[] markArray=new int[solution.length];//同一类的就为同一个标记
        int curMark=1;//0的标记为0
        List<List<Integer>>vehicleList=splitVehicle(solution);
        for (int i = 0; i < vehicleList.size(); i++) {
            List<Integer> vehicle=vehicleList.get(i);
            if(vehicle.size()>=K){
                List<List<Integer>>clusters=kMeans(points,distanceMatrix,vehicle);
                for (int j = 0; j < clusters.size(); j++) {
                    List<Integer>clusterJ=clusters.get(j);
                    if(clusterJ.size()>=2){
                        for (int k = 0; k < clusterJ.size(); k++) {
                            markArray[clusterJ.get(k)]=curMark;
                        }
                        curMark++;
                    }
                }
            }
        }
        if (curMark>=3){
            // 获取每个mark的索引
            Map<Integer, List<Integer>> markIndicesMap = getMarkIndices(solution, markArray);
            int i=0,tryTime=300;
            while(i<tryTime){
                int randMark=r.nextInt(curMark-2)+1;//0的不换
                 swapRandomElements(solution,markIndicesMap,randMark);
                 solution[solution.length-1]=cucalateTotalLength(solution,distanceMatrix);
                 newSolutions.add(Arrays.copyOf(solution, solution.length));
                 i++;
            }

        }
        sortSolutionsByTotalLength(newSolutions);
        double maxValue = bestValue * 1.015;
        int count = 0;
        int lengthIndex = newSolutions.get(0).length - 1;

        while (count < newSolutions.size() && newSolutions.get(count)[lengthIndex] <=maxValue) {
            count++;
        }

        newSolutions.subList(count, newSolutions.size()).clear();
        return newSolutions;


    }
    // 获取每个mark对应的索引列表
    public static Map<Integer, List<Integer>> getMarkIndices(double[] solution, int[] markArray) {
        Map<Integer, List<Integer>> markIndicesMap = new HashMap<>();

        for (int i = 0; i < markArray.length; i++) {
            int mark = markArray[i];
            markIndicesMap.putIfAbsent(mark, new ArrayList<>());
            markIndicesMap.get(mark).add(i);
        }

        return markIndicesMap;
    }
    public static void swapRandomElements(double[] solution, Map<Integer, List<Integer>> markIndicesMap, int targetMark) {
        List<Integer> indices = markIndicesMap.get(targetMark);

        Random random = new Random();
        // 随机选择两个不同的索引
        int index1 = random.nextInt(indices.size());
        int index2;
        do {
            index2 = random.nextInt(indices.size());
        } while (index1 == index2);

        // 交换两个元素
        int swapIndex1 = indices.get(index1);
        int swapIndex2 = indices.get(index2);

        // 进行交换
        double temp = solution[swapIndex1];
        solution[swapIndex1] = solution[swapIndex2];
        solution[swapIndex2] = temp;

    }

    //对车辆里的客户进行聚类，分为3类
    public static List<List<Integer>> kMeans(double[][] points, double[][] distanceMatrix, List<Integer> vehicle) {
        int n = vehicle.size();
        double[][] centroids = new double[K][2];
        List<List<Integer>> clusters = new ArrayList<>();

        // 随机选择初始质心
        for (int i = 0; i < K; i++) {
            centroids[i] = new double[]{points[vehicle.get(i)][0], points[vehicle.get(i)][1]}; // 选择前两列作为初始质心
            clusters.add(new ArrayList<>());
        }

        boolean changed;
        int iteration = 0;

        do {
            // 清空每个簇
            for (List<Integer> cluster : clusters) {
                cluster.clear();
            }

            // 将每个点分配给最近的质心
            for (int i = 0; i < n; i++) {
                int closestCentroid = findClosestCentroid(vehicle.get(i), centroids,points);
                clusters.get(closestCentroid).add(vehicle.get(i));
            }

            changed = updateCentroids(clusters, centroids, points);
            iteration++;
        } while (changed && iteration < MAX_ITERATIONS); // 如果质心发生变化且未达到最大迭代次数，继续迭代

        return clusters;
    }

    private static int findClosestCentroid(int pointIndex, double[][] centroids,double[][]points) {
        int closest = 0;
        double minDistance = Double.MAX_VALUE;

        for (int i = 0; i < centroids.length; i++) {
            double distance = Math.pow(points[pointIndex][0] -centroids[i][0], 2) + Math.pow(points[pointIndex][1] -centroids[i][1], 2);
            if (distance < minDistance) {
                minDistance = distance;
                closest = i;
            }
        }

        return closest;
    }

    private static boolean updateCentroids(List<List<Integer>> clusters, double[][] centroids, double[][] points) {
        boolean changed = false;

        for (int i = 0; i < K; i++) {
            if (clusters.get(i).isEmpty()) continue;
            double[] newCentroid = calculateCentroid(clusters.get(i), points);
            if (!equals(centroids[i], newCentroid)) {
                centroids[i] = newCentroid;
                changed = true;
            }
        }

        return changed;
    }

    private static double[] calculateCentroid(List<Integer> cluster, double[][] points) {
        double sumX = 0;
        double sumY = 0;

        for (int cust : cluster) {
            sumX += points[cust][0]; // 只取坐标X
            sumY += points[cust][1]; // 只取坐标Y
        }

        return new double[]{sumX / cluster.size(), sumY / cluster.size()};
    }

    private static boolean equals(double[] a, double[] b) {
        return a[0] == b[0] && a[1] == b[1];
    }



}