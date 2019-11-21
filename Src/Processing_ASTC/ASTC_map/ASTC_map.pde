/*Accumulative SpatioTemporal Contraction map //<>//
This code plots the ASTC map by importing xyt triplets from a csv file.
This csv file may be created from the output of export_tree_centroids in
the AR-DIC Toolbox in Matlab.
Modify the file name in loadTable to open other csv files.

All paths are plotted in gray except those specified in the vector "nums"
which are modified for specific colors.

Map controls:
Left click and drag to rotate
Right click and drag or roll center mouse button to zoom
Click and drag with center mouse button to pan

Adaptive Reference Digital Image Correlation v 1.0 2018
Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
*/

//Initialize and setup:
import peasy.*;// include PeasyCam library
PeasyCam cam; //initialize variables
int[] nums={10, 11, 21, 19};//specify which path numbers to color
float[][] x;// array to contain x coordinates
float[][] y;//array to contain y coordinates
float[][] t;//array to contain time coordinates
Table table;// table for storing values from csv file
int numpath;
int cols;
int subindex;
int R;
int G;
int B;
float t_start=0;//set time minimum dimension
float t_end=1000;//set time maximum dimension
float t_temp;
float t_temp2;

PImage bg;
void setup() {
  size(600, 600, P3D);//window size
  bg=loadImage("background.png");//load background image
  cam = new PeasyCam(this, 20, 20, 100, 150); //create camera
  cam.rotateX(-3.14/7);//set camera view
  table = loadTable("ASTC_example_data.csv"); //load table from csv file. xyt triplets
  numpath=table.getRowCount(); // each row is a tree path
  cols=table.getColumnCount();
  x = new float[int(cols/3)][numpath];
  y = new float[int(cols/3)][numpath];
  t = new float[int(cols/3)][numpath];
  for (int row=0; row<numpath; row++) { //translate table into x,y,t arrays
    TableRow  r=table.getRow(row);
    subindex=0;
    for (int i=0; i<r.getColumnCount(); i=i+3) {
      x[subindex][row]=r.getFloat(i);
      y[subindex][row]=r.getFloat(i+1);
      t[subindex][row]=r.getFloat(i+2);
      subindex++;
    }
  }
  R=250;
  B=0;
  G=0;
  sphereDetail(30);
}//setup

void draw() {
  background(bg);// set image as background
  lights();
  spotLight(255, 255, 255, width/8, height/8, 400, 0, 0, -1, PI/4, 2);
  fill(0, 0, 0);
  noStroke();
  specular(255, 255, 255);
  shininess(100);
  rect(0, 0, width/8, height/8);
  noFill();
  for (int j=0; j<numpath; j++) { //draw paths
    R=73; //set nodes not origins to gray
    G=73;
    B=73;
    stroke(R, G, B);//set specified stroke color
    strokeWeight(2);
    beginShape();
    for (int i=0; i<cols/3; i++) {
      t_temp=t[i][j];
      if ((t_temp>t_start)&(t_temp<t_end)) {
        curveVertex(x[i][j], y[i][j], t_temp);
      }
    }
    endShape();
    strokeWeight(8);

    for (int ii=0; ii<cols/3; ii++) {//draw spheres at nodes
      t_temp2=t[ii][j];
      if ((t_temp2>t_start)&(t_temp2<t_end)) {
        fill(R, G, B);
        noStroke();
        pushMatrix();
        translate(x[ii][j], y[ii][j], t_temp2);
        sphere(1.5);
        popMatrix();
      }
      noFill();
      stroke(R, G, B);
    }
  }

  for (int j : nums) { //draw paths with distinct color
    if (j==10) {//0'''
      R=244;
      G=195;
      B=137;
    } else if (j==11) {//O''' second half
      R=244;
      G=195;
      B=137;
    } else if (j==19) {//O'
      R=121;
      G=20;
      B=54;
    } else if (j==21) {//O''
      R=2;
      G=38;
      B=60;
    }       
    stroke(R, G, B);//set specified stroke colors
    strokeWeight(2);
    beginShape();
    for (int i=0; i<cols/3; i++) {
      t_temp=t[i][j];
      if ((t_temp>t_start)&(t_temp<t_end)) {
        curveVertex(x[i][j], y[i][j], t_temp);
      }
    }
    endShape();
    strokeWeight(8);
    for (int ii=0; ii<cols/3; ii++) {//draw nodes as spheres with specific colors
      t_temp2=t[ii][j];
      if ((t_temp2>t_start)&(t_temp2<t_end)) {
        fill(R, G, B);
        noStroke();
        pushMatrix();
        translate(x[ii][j], y[ii][j], t_temp2);//for 2D appearance: point(x[ii][j], y[ii][j], t_temp2);
        sphere(1.5);
        popMatrix();
      }
      noFill();
      stroke(R, G, B);
    }
  }

  drawAxis();
}//draw

void drawAxis() {
  stroke(255, 255, 255);
  strokeWeight(4);
  textSize(20);
  line(0, 0, 0, 56, 0, 0);//x axis
  line(0, 0, 0, 0, 41, 0);//yaxis
  line(0, 0, 0, 0, 0, 160);//zaxis
  fill(255, 255, 255);
  text('x', 56, 0, 0);
  text('y', -12, 41, 0);
  rotateX(-3.14/2);
  fill(255, 255, 255);
  translate(-10, -160, 0);
  text('t', 0, 0, 0);
}