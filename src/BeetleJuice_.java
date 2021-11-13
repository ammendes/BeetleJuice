import ij.IJ;
import ij.plugin.PlugIn;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import java.io.*;
import java.util.Random;

import static java.lang.Math.*;

public class BeetleJuice_ implements PlugIn {

    //---- Define variables ----
    int windowWidth = 2500; // Width of the FOV (in pixels)
    int windowHeight = 2500; // Height of the FOV (in pixels)
    int nFrames = 20000; // Total number of frames to simulate
    float blinksPerFrame = 0.041F; // Empirical; 820 localisations in a ~200 nm diameter Gag VLP with STORM (55% DoL).
    int nBlinks = (int) (nFrames*blinksPerFrame);
    int framesPerBlink = nFrames/nBlinks;
    String[] localisationTable = new String[(int) nBlinks+1]; // +1 because first index is the table headers
    static int id = 1; // The ID of each blinking event (starts at 1 and will be incremented dynamically)
    final int initialX = windowWidth/2; // particle is centered in X
    final int initialY = windowHeight/2; // particle is centered in Y
    final int maxRadius = 75; // Radius of the circle representing the vesicle (in nanometers?)

    // Empirical parameters taken from a virus I imaged with STORM
    float sigmaMean = 110;
    float sigmaVar = 565;
    float sigmaMin = 50;
    float sigmaMax = 155;
    float intensityMean = 660;
    float intensityStdDev = 450;
    float intensityMin = 250;
    float intensityMax = 2740;
    float offsetMean = 470;
    float offsetVar = 1330;
    float offsetMin = 410;
    float offsetMax = 630;
    float bkgStdMean = 30;
    float bkgStdVar = 16;
    float bkgStdMin = 20;
    float bkgStdMax = 50;
    float chi2Mean = 100;
    float chi2Var = 400;
    float chi2Min = 50;
    float chi2Max = 200;
    float uncertaintyMean = 20;
    float uncertaintyVar = 25;
    float uncertaintyMin = 10;
    float uncertaintyMax = 30;


    @Override
    public void run(String s) {
        //---- Get blinks for each frame ----
        for (int frame = 1; frame <= nFrames; frame = frame + framesPerBlink) {
            GetBlinks blinks = new GetBlinks();
            blinks.setup(localisationTable, id, frame, initialX, initialY, maxRadius, sigmaMean, sigmaVar, sigmaMin,
                    sigmaMax, intensityMean, intensityStdDev, intensityMin, intensityMax, offsetMean, offsetVar,
                    offsetMin,offsetMax, bkgStdMean, bkgStdVar, bkgStdMin, bkgStdMax, chi2Mean, chi2Var, chi2Min,
                    chi2Max, uncertaintyMean, uncertaintyVar, uncertaintyMin, uncertaintyMax);
            blinks.start();
            id += 1;
        }

        // Add headers to the final table and write to file (.csv)
        localisationTable[0] = "id, frame, x [nm], y [nm], sigma [nm], intensity [photon], offset [photon], bkgstd [photon], chi2, uncertainty [nm]";

        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter("/Users/ammendes/Desktop/simulation.csv"));

            for(int i=0; i<localisationTable.length;i++){
                bw.write(localisationTable[i]);
                bw.newLine();
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        IJ.log("Done!");
    }

}
    
class GetBlinks extends Thread {
    String[] localisationTable;
    int id, frame, initialX, initialY;
    float maxRadius, sigma, intensity, offset, bkgStd, chi2, uncertainty;
    float[] xy;

    public void setup(String[] localisationTable, int id, int frame, int initialX, int initialY, float maxRadius,
                      float sigmaMean, float sigmaVar, float sigmaMin, float sigmaMax, float intensityMean,
                      float intensityStdDev, float intensityMin, float intensityMax, float offsetMean, float offsetVar,
                      float offsetMin, float offsetMax, float bkgStdMean, float bkgStdVar, float bkgStdMin,
                      float bkgStdMax, float chi2Mean, float chi2Var, float chi2Min, float chi2Max,
                      float uncertaintyMean, float uncertaintyVar, float uncertaintyMin, float uncertaintyMax)
    {

        this.id = id;
        this.frame = frame;
        this.initialX = initialX;
        this.initialY = initialY;
        this.maxRadius = maxRadius;
        this.localisationTable = localisationTable;
        xy = getXY(initialX, initialY, maxRadius);
        sigma = getTruncatedNormal(sigmaMean, sigmaVar, sigmaMin, sigmaMax);
        intensity = getLogNormal(intensityStdDev, intensityMean, intensityMin, intensityMax);
        offset = getTruncatedNormal(offsetMean, offsetVar, offsetMin, offsetMax);
        bkgStd = getTruncatedNormal(bkgStdMean, bkgStdVar, bkgStdMin, bkgStdMax);
        chi2 = getTruncatedNormal(chi2Mean, chi2Var, chi2Min, chi2Max);
        uncertainty = getTruncatedNormal(uncertaintyMean, uncertaintyVar, uncertaintyMin, uncertaintyMax);
    }

    @Override
    public void run() {
        // Store localisations in the final table
        localisationTable[id] = id + "," + frame + "," + xy[0] + "," + xy[1] + "," + sigma + "," + intensity + "," + offset +
                "," + bkgStd + "," + chi2 + "," + uncertainty;
    }

    //---- USER METHODS ----
    // Get XY position (uniformly distributed points in a circle of defined radius)
    public float[] getXY(float initialX, float initialY, float maxRadius) {
        float[] positions = new float[2];
        Random rand = new Random();

        for(int i = 1; i<=1500; i++) {
            float stepMagnitude;
            float theta = (float) (rand.nextFloat() * 2 * PI);
            float u = rand.nextFloat() + rand.nextFloat();
            if (u > 1) {
                stepMagnitude = 2 - u;
            } else {
                stepMagnitude = u;
            }
            positions[0] = (float) (maxRadius*stepMagnitude * cos(theta) + initialX);
            positions[1] = (float) (maxRadius*stepMagnitude * sin(theta) + initialY);
        }
        return positions;
    }




    // Get Y position (based on random motion)
    public float getY(float initialY) {
        Random rand = new Random();
        float theta = (float) (2 * PI * rand.nextFloat());
        float amplitude = (float) sqrt(pow(maxRadius, 2) * rand.nextFloat());

        return (float) (initialY+amplitude * sin(theta));
    }

        // Points are too centered
        /*
        Random rand = new Random();
        float a = rand.nextFloat();
        float b = rand.nextFloat();
        float nextPosition;

        if (a<b) {
            nextPosition = (float) (initialY+b*maxRadius*sin(2*PI*(a/b)));
        }else{
            nextPosition = (float) (initialY+b*maxRadius*sin(2*PI*(b/a)));
        }
        return nextPosition;
        */

        // More uniform than the previous one
        /*
        Random rand = new Random();
        float amplitude = (float) (maxRadius * sqrt(rand.nextFloat())); // Amplitude
        float alpha = (float) (rand.nextFloat()*2*PI); // Angle (in radians)

        return (float) round(amplitude*sin(alpha)+initialY);
        */


    float getTruncatedNormal(float mean, float variance, float min, float max) {
        Random rand = new Random();
        double number = mean + rand.nextGaussian() * variance;
        while (number < min || number > max){
            number = mean + rand.nextGaussian() * variance;
        }
        return (float) number;
    }

    // Generate random numbers from a Log-Normal distribution
    float getLogNormal(float sigma, float mean, float min, float max){
        LogNormalDistribution logNormal = new LogNormalDistribution(sigma, mean);
        double number = logNormal.sample();
        while (number > max || number < min){
            number = logNormal.sample();
        }
        return (float) number;
    }
}
