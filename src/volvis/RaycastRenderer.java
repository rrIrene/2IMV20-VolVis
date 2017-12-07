/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import static java.lang.Math.sqrt;
import java.util.Arrays;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    private void mip(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        double diagonal = VectorMath.length(new double[]{volume.getDimX(),volume.getDimY(),volume.getDimZ()});
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {          
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + viewVec[0] * volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + viewVec[1] * volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + viewVec[2] * volumeCenter[2];
                
                int maxVoxel = 0;
                
                //pixelCoord is the furthest point on the ray between this pixel through the volume data that is still within the bounding box.
                //iterate through the ray by following the unit vector towards the viewing plane in order to find the maximum voxel.
                for (int step = 0; step < Math.floor(diagonal); step++) {
                    try {
                       int val = getVoxel(pixelCoord);
                       if (maxVoxel < val)  maxVoxel = val;
                       pixelCoord[0] -= viewVec[0];
                       pixelCoord[1] -= viewVec[1];
                       pixelCoord[2] -= viewVec[2];                        
                    } catch (Exception ex) {
                        System.out.println("Exception at step: " + step + ": " + Arrays.toString(pixelCoord));
                    }
                }
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxVoxel/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxVoxel > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(maxVoxel);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }
    
    private void compositing(double[] viewMatrix) {
        //Use raytracing similar to MIP
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        double diagonal = VectorMath.length(new double[]{volume.getDimX(),volume.getDimY(),volume.getDimZ()});
        
        // sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();
        //stepVector has the same direction as viewVec, but has length of stepLength.
        double[] stepVector = new double[3];
        double stepLength = 1;
        if (diagonal > compositingStep) stepLength = Math.floor(diagonal / compositingStep);
        VectorMath.setVector(stepVector, viewVec[0] * stepLength, viewVec[1] * stepLength, viewVec[2] * stepLength);
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {          
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + viewVec[0] * volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + viewVec[1] * volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + viewVec[2] * volumeCenter[2];
                   
                double voxelval = Interpolate(pixelCoord);
                int val = (int) voxelval;  
                //int val = getVoxel(pixelCoord);
                voxelColor = tFunc.getColor(val);
                double C_r = voxelColor.r;
                double C_g = voxelColor.g;
                double C_b = voxelColor.b;
                double C_a = voxelColor.a;
                
                //pixelCoord is the furthest point on the ray between this pixel through the volume data that is still within the bounding box.
                //iterate through the ray by following the unit vector towards the viewing plane in order to calculate the color and opacity of this pixel.
                for (double step = stepLength; step < Math.floor(diagonal); step+=stepLength) {
                    try {
                       voxelval = Interpolate(pixelCoord);
                       val = (int) voxelval;
                       //val = getVoxel(pixelCoord);
                       voxelColor = tFunc.getColor(val);
                       C_r = (1 - voxelColor.a) * C_r + voxelColor.a * voxelColor.r;
                       C_g = (1 - voxelColor.a) * C_g + voxelColor.a * voxelColor.g;
                       C_b = (1 - voxelColor.a) * C_b + voxelColor.a * voxelColor.b;
                       C_a = (1 - voxelColor.a) * C_a + voxelColor.a;
                       pixelCoord[0] -= stepVector[0];
                       pixelCoord[1] -= stepVector[1];
                       pixelCoord[2] -= stepVector[2];                        
                    } catch (Exception ex) {
                        System.out.println("Exception at step: " + step + ": " + Arrays.toString(pixelCoord));
                    }
                }
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = C_a <= 1.0 ? (int) Math.floor(C_a * 255) : 255;
                int c_red = C_r <= 1.0 ? (int) Math.floor(C_r * 255) : 255;
                int c_green = C_g <= 1.0 ? (int) Math.floor(C_g * 255) : 255;
                int c_blue = C_b <= 1.0 ? (int) Math.floor(C_b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }

    public void setCompositingStep(int step) {
        compositingStep = step;
    }

    private void tf2d(double[] viewMatrix) {
        //Use raytracing similar to MIP
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        double diagonal = VectorMath.length(new double[]{volume.getDimX(),volume.getDimY(),volume.getDimZ()});
        
        // sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();
        //stepVector has the same direction as viewVec, but has length of stepLength.
        double[] stepVector = new double[3];
        double stepLength = 1;
        if (diagonal > compositingStep) stepLength = Math.floor(diagonal / compositingStep);
        VectorMath.setVector(stepVector, viewVec[0] * stepLength, viewVec[1] * stepLength, viewVec[2] * stepLength);
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {          
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + viewVec[0] * volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + viewVec[1] * volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + viewVec[2] * volumeCenter[2];
               
                //Calculate the color & opacity using tf2d mapping
                voxelColor = getTF2dColor(pixelCoord, viewVec);
                double C_a = 1 - voxelColor.a;
                double C_r = voxelColor.r;
                double C_g = voxelColor.g;
                double C_b = voxelColor.b;

                //pixelCoord is the furthest point on the ray between this pixel through the volume data that is still within the bounding box.
                //iterate through the ray by following the unit vector towards the viewing plane in order to calculate the color and opacity of this pixel.
                for (double step = stepLength; step < Math.floor(diagonal); step+=stepLength) {
                    voxelColor = getTF2dColor(pixelCoord, viewVec);
                    C_a = (1 - voxelColor.a) * C_a;
                    if (volumeShading) {
                        //Only reconsider voxelColor when we're applying volume shading. Otherwise, just use the color selected by triangleWidget.
                        C_r = (1 - voxelColor.a) * C_r + voxelColor.a * voxelColor.r;
                        C_g = (1 - voxelColor.a) * C_g + voxelColor.a * voxelColor.g;
                        C_b = (1 - voxelColor.a) * C_b + voxelColor.a * voxelColor.b;
                    }
                    pixelCoord[0] -= stepVector[0];
                    pixelCoord[1] -= stepVector[1];
                    pixelCoord[2] -= stepVector[2];  
                }           
                
                C_a = 1 - C_a;
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = C_a <= 1.0 ? (int) Math.floor(C_a * 255) : 255;
                int c_red = C_r <= 1.0 ? (int) Math.floor(C_r * 255) : 255;
                int c_green = C_g <= 1.0 ? (int) Math.floor(C_g * 255) : 255;
                int c_blue = C_b <= 1.0 ? (int) Math.floor(C_b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }

    private TFColor getTF2dColor(double[] coord, double[] V) {
        //Get important values from tf2d editor
        int f_v = tfEditor2D.triangleWidget.baseIntensity;
        double r = tfEditor2D.triangleWidget.radius;
        
        int f_x = getVoxel(coord);
        VoxelGradient grad = getGradient(coord);
        double delta_fx = grad.mag;
        
        //check if the magnitude is within range.
        boolean withinRange = delta_fx >= tfEditor2D.triangleWidget.range1 && delta_fx <= tfEditor2D.triangleWidget.range2;
        
        double factor = 0;
        
        //Get the opacity value using Levoy's formula.
        if (delta_fx == 0.0 && f_v == f_x) {
            factor = 1;
        } else if (delta_fx > 0.0 && f_x-(r*delta_fx) <= f_v && f_v <= f_x+(r*delta_fx) && withinRange) {
            factor = 1 - Math.abs((f_v-f_x)/delta_fx) / r;
        }
        
        TFColor color_v = new TFColor();
        color_v.a = tfEditor2D.triangleWidget.color.a * factor;
        color_v.r = tfEditor2D.triangleWidget.color.r;
        color_v.g = tfEditor2D.triangleWidget.color.g;
        color_v.b = tfEditor2D.triangleWidget.color.b; 
        
        //implement Phong shading to get the right colors, if applicable
        if (volumeShading) {
            //Kniss: "The normalized gradient is often used as the normal for surface-based volume shading"
            double[] N = new double[3];
            VectorMath.setVector(N, grad.x/grad.mag, grad.y/grad.mag, grad.z/grad.mag);
            //Calculate the dot products
            double diffProduct = VectorMath.dotproduct(V, N);
            if (diffProduct > 0.0) {
                double specProduct = VectorMath.dotproduct(N, V);   //H = L+V/|L+V| = V, because V = viewVec which is already normalized.
                if (specProduct > 0.0) {
                    color_v.r = tfEditor2D.ambientColor.r * tfEditor2D.k_ambient + tfEditor2D.triangleWidget.color.r * tfEditor2D.k_diff * diffProduct + tfEditor2D.k_spec * Math.pow(specProduct, tfEditor2D.phong_alpha);
                    color_v.g = tfEditor2D.ambientColor.g * tfEditor2D.k_ambient + tfEditor2D.triangleWidget.color.g * tfEditor2D.k_diff * diffProduct + tfEditor2D.k_spec * Math.pow(specProduct, tfEditor2D.phong_alpha);
                    color_v.b = tfEditor2D.ambientColor.b * tfEditor2D.k_ambient + tfEditor2D.triangleWidget.color.b * tfEditor2D.k_diff * diffProduct + tfEditor2D.k_spec * Math.pow(specProduct, tfEditor2D.phong_alpha);
                }
            }
        } 
             
        return color_v;
    }
    
    public enum RaycastRenderType {
        SLICER, MIP, COMPOSITING, TF2D
    }
    
    private int compositingStep = 1;
    private RaycastRenderType type = RaycastRenderType.SLICER;
    private boolean volumeShading = false;
    
    public void setRenderType(RaycastRenderType type) {
        this.type = type;
    }
    
    public void toggleVolumeShading() {
	this.volumeShading = !this.volumeShading;
    }
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
                
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }
    
    VoxelGradient getGradient(double[] coord) {
        if (coord[0] < 0 || coord[0] >= gradients.getDimX() || coord[1] < 0 || coord[1] >= gradients.getDimY()
                || coord[2] < 0 || coord[2] >= gradients.getDimZ()) {
            return new VoxelGradient();
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        return gradients.getGradient(x, y, z);
    }

    void slicer(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0]-0.5;
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1]-0.5;
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2]-0.5;
                
                double voxelval = this.Interpolate(pixelCoord);
                int val = (int) voxelval;
//                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                //voxelColor.r = val/max;
                //voxelColor.g = voxelColor.r;
                //voxelColor.b = voxelColor.r;
                //voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
        

    }

    private double Interpolate(double[] pixelCoord) {        
         
        if (pixelCoord[0] < 0 || pixelCoord[0] > volume.getDimX()-1 || pixelCoord[1] < 0 || pixelCoord[1] > volume.getDimY()-1
                || pixelCoord[2] < 0 || pixelCoord[2] > volume.getDimZ()-1) {
            return 0;
        }
        else{

         double[] C000 = new double[] { Math.floor(pixelCoord[0]), 
                                Math.floor(pixelCoord[1]),
                                Math.floor(pixelCoord[2]) };
         double[] C100 = new double[] { Math.ceil(pixelCoord[0]), 
                                Math.floor(pixelCoord[1]),
                                Math.floor(pixelCoord[2]) };
         double[] C001 = new double[] {Math.floor(pixelCoord[0]), 
                                 Math.floor(pixelCoord[1]),
                                Math.ceil(pixelCoord[2]) };
         double[] C101 = new double[] { Math.ceil(pixelCoord[0]), 
                                Math.floor(pixelCoord[1]),
                                Math.ceil(pixelCoord[2]) };
         double[] C010 = new double[] {Math.floor(pixelCoord[0]), 
                                Math.ceil(pixelCoord[1]),
                                Math.floor(pixelCoord[2]) };
       
         double[] C110 = new double[] { Math.ceil(pixelCoord[0]), 
                                Math.ceil(pixelCoord[1]),
                                Math.floor(pixelCoord[2]) };
       
         double[] C011 = new double[] {Math.floor(pixelCoord[0]), 
                                Math.ceil(pixelCoord[1]),
                                Math.ceil(pixelCoord[2]) };
         double[] C111 = new double[] {Math.ceil(pixelCoord[0]), 
                                Math.ceil(pixelCoord[1]),
                                Math.ceil(pixelCoord[2]) };
       
        
        double alpha = (pixelCoord[0] -(int) Math.floor(pixelCoord[0]))/((int) Math.ceil(pixelCoord[0]) - (int) Math.floor(pixelCoord[0]));
        double beta = (pixelCoord[1] -(int) Math.floor(pixelCoord[1]))/((int) Math.ceil(pixelCoord[1]) - (int) Math.floor(pixelCoord[1]));
        double gama = (pixelCoord[2] -(int) Math.floor(pixelCoord[2]))/((int) Math.ceil(pixelCoord[2]) - (int) Math.floor(pixelCoord[2]));
        
        
        if(Double.isNaN(alpha)) alpha = 0.0;
        if(Double.isNaN(beta))  beta = 0.0;
        if(Double.isNaN(gama))  gama = 0.0;
        
        double C00 = (1-alpha)*getVoxel(C000) + alpha*getVoxel(C100); 
        double C01 = (1-alpha)*getVoxel(C001) + alpha*getVoxel(C101);
        double C10 = (1-alpha)*getVoxel(C010) + alpha*getVoxel(C110);
        double C11 = (1-alpha)*getVoxel(C011) + alpha*getVoxel(C111);
        
        double C0 = (1-beta)*C00 + beta*C10;
        double C1 = (1-beta)*C01 + beta*C11;
        double C = (1-gama)*C0 + gama*C1;
           
        
        return C;
        }
    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        
        switch(type) {
            case SLICER:
                slicer(viewMatrix); 
                break;
            case MIP:
                mip(viewMatrix);
                break;
            case COMPOSITING:
                compositing(viewMatrix);
                break;
            case TF2D:
                tf2d(viewMatrix);
                break;
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
