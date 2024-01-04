#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <math.h>
#include <set>
#include <ctime>
#include <cstdlib>
#include <random>
std::random_device rd;
std::mt19937 gen(rd());

class Pixel{
private:
    
    std::uniform_real_distribution<double> distribution;
    double diff(double val, double pert){
        double to_return = val-pert;
        if(to_return<0){
            return 0;
        }
        return to_return;
    }
    double not_diff(double val, double pert){
        double to_return = val+pert;
        if(to_return>255){
            return 255;
        }
        return to_return;
    }

public:
    double r;
    double g;
    double b;
    Pixel(double r1, double g1, double b1){
        this->r = r1;
        this->g = g1;
        this->b = b1;
    }

    static double distance(Pixel pix1, Pixel pix2){
        return fabs(pix1.r-pix2.r)+fabs(pix1.b-pix2.b)+fabs(pix1.g-pix2.g);

    }

    Pixel roundPix(){
        double r = (double)round(this->r);
        double g = (double)round(this->g);
        double b = (double)round(this->b);
        return Pixel(r,g,b);
    }

    std::pair<Pixel, Pixel> perturb(){
        double rand1 = distribution(gen);
        double rand2 = distribution(gen);
        double rand3 = distribution(gen);

        double r1 = this->r;
        double b1 = this->b;
        double g1 = this->g;

        double r2 = this->r;
        double b2 = this->b;
        double g2 = this->g;

        Pixel p1 = Pixel(diff(r1, rand1),diff(g1, rand2), diff(b1, rand3));
        Pixel p2 = Pixel(not_diff(r1, rand1),not_diff(g1, rand2), not_diff(b1, rand3));
        return std::make_pair(p1, p2);
    }

    bool operator<(const Pixel& other) const {
        if (r != other.r) {
            return r < other.r;
        } else if (g != other.g) {
            return g < other.g;
        } else {
            return b < other.b;
        }
    }
};

class Quantizer{
    static Pixel avg_pix(std::vector<Pixel> to_avg){
        double r = 0;
        double g = 0;
        double b = 0;
        for(int i =0; i<to_avg.size(); i++){
            r+=to_avg[i].r;
            g+=to_avg[i].g;
            b+=to_avg[i].b;
        }
        return Pixel(r/to_avg.size(), g/to_avg.size(), b/to_avg.size());
    }
public:
    int uniquePixelsCount;
    std::vector<Pixel> pixels;
    std::vector<Pixel> codebookMain;
    std::vector<Pixel> outputMain;
    int colors;
    Quantizer(std::vector<Pixel> pixels1, int colors1){
        this->pixels = pixels1;
        std::set<Pixel> uniquePixels(pixels1.begin(), pixels1.end());
        this->uniquePixelsCount = uniquePixels.size();
        this->colors=colors1;
        this->codebookMain = quantize();
        this->outputMain = to_output();
    }
    
    std::vector<Pixel> quantize(){
        std::vector<Pixel> codebook;
        Pixel mean = avg_pix(this->pixels);
        codebook.push_back(mean);
        int sizeTarget = std::pow(2, this->colors);
        if(sizeTarget>this->uniquePixelsCount){
            printf("Quantization is useless, this image has %d unique colors, you demand %d colors", uniquePixelsCount, sizeTarget);
            return this->pixels;
        }
        //for(int i = 0; i<this->colors; i++){
        while(sizeTarget>codebook.size()){
            printf("codebook size: %d\n", codebook.size());
            int beforePerturb = codebook.size();
            printf("target size: %d\n", sizeTarget);
            std::vector<Pixel> new_codebook;
            for (int j = 0; j<codebook.size(); j++){
                auto pixels = codebook[j].perturb();
                new_codebook.push_back(pixels.first);
                new_codebook.push_back(pixels.second);
            }
            codebook = new_codebook;
            long double prev_dist=0;
            while(true){
                std::vector<std::vector<Pixel> > clusters;
                for (int j=0; j<codebook.size(); j++){
                    std::vector<Pixel> p;
                    clusters.push_back(p);
                }

                for(int j=0; j<this->pixels.size(); j++){
                    int minDist = INT32_MAX;
                    int index = INT32_MAX;
                    int k;
                    for(k=0; k<codebook.size(); k++){
                        int dist = Pixel::distance(codebook[k], this->pixels[j]);
                        if(dist<minDist){
                            minDist = dist;
                            index = k;
                        }
                    }
                    clusters[index].push_back(this->pixels[j]);
                }
                
                long double distortion =0;
                for(int j=0; j<clusters.size(); j++){
                    for(int k=0; k<clusters[j].size(); k++){
                        distortion+=Pixel::distance(codebook[j], clusters[j][k]);
                    }
                }
                //printf("distorion: %d\n", distortion);
                if(distortion==0){
                    //printf("distorion: %d\n", distortion);
                    break;
                }else if(fabs(distortion-prev_dist)/distortion < 0.001){
                    //printf("distorion: %d\n", distortion);
                    break;
                }
                prev_dist = distortion;
                std::vector<Pixel> new_centroids;
                for(int j=0; j<codebook.size(); j++){
                    if(!clusters[j].empty()){
                    Pixel p = avg_pix(clusters[j]);
                    new_centroids.push_back(p);
                    }
                }
                codebook = new_centroids;
            }
            if(codebook.size()==beforePerturb){
                //zabezpieczenie przed wpadnięciem w minima
                return codebook;
            }
        }
        //}
        return codebook;
    }

    std::vector<Pixel> to_output(){
        std::vector<Pixel> output;
        for(int i=0; i<this->pixels.size(); i++){
            int minDist=INT32_MAX;
            int index = -1;
            for(int j=0; j<this->codebookMain.size(); j++){
                int dist = Pixel::distance(this->pixels[i], this->codebookMain[j]);
                if(dist<minDist){
                    index = j;
                    minDist=dist;
                }
            }
            output.push_back(this->codebookMain[index].roundPix());
        }
        double mse=0;
        double snr=0;
        for(int i=0; i<this->pixels.size(); i++){
            mse+=std::pow(Pixel::distance(output[i], this->pixels[i]), 2.0);
            snr+=std::pow(this->pixels[i].r, 2)+std::pow(this->pixels[i].g, 2)+std::pow(this->pixels[i].b, 2);
        }
        mse/=static_cast<double>(this->codebookMain.size());
        snr/=static_cast<double>(this->codebookMain.size());
        snr/=mse;
        printf("MSE = %f\n", mse);
        printf("SNR = %f\n", snr);
        return output;
    }

};

class Image{
public:
    std::string sourceFile;
    std::string destFile;
    std::vector<uint8_t> inputBytes;
    std::vector<Pixel> output;

    Image(std::string source, std::string destination, int colors){
        this->sourceFile = source;
        this->destFile = destination;
        this->inputBytes = tgaBytes();
        this->output = Quantizer(tgaPixels(), colors).outputMain;
        writeTga();
    }

    std::vector<uint8_t> tgaBytes(){
        std::ifstream stream(this->sourceFile, std::ios::in | std::ios::binary);
        std::vector<uint8_t> contents((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
        return contents;
    }

    std::vector<Pixel> tgaPixels(){
        std::vector<Pixel> pixels;
        int b12 = this->inputBytes[12];
        int b13 = this->inputBytes[13];
        int b14 = this->inputBytes[14];
        int b15 = this->inputBytes[15];
        int width = b13*256+b12;
        int height = b15*256+b14;
        int startPixels = 18;
        int endPixels = width*height*3+18;
        for(int i=startPixels; i<endPixels; i+=3){
            double r = this->inputBytes[i];
            double g = this->inputBytes[i + 1];
            double b = this->inputBytes[i + 2];
            Pixel p = Pixel(r,g,b);
            pixels.push_back(p);
        }
        return pixels;
    }
    void writeTga(){
        int byteCounter = 0;
        std::ofstream outFile(this->destFile, std::ios::binary);
        if (!outFile.is_open()) {
            std::cerr << "Couldn't open file for writing: " << this->destFile << std::endl;
            return;
        }
        for(int i=0; i<18;i++){
            outFile.put(this->inputBytes[i]);
            byteCounter++;
        }
        for(int i=0; i<output.size();i++){
            unsigned char r = output[i].r;
            unsigned char g = output[i].g;
            unsigned char b = output[i].b;
            outFile.put(r);
            outFile.put(g);
            outFile.put(b);
            byteCounter+=3;
        }
        for(int i=byteCounter; i<inputBytes.size(); i++){
            outFile.put(this->inputBytes[i]);
        }
    }
};

int main(int argc, char* argv[]){
    if(argc!=4){
        printf("Usage: <source file> <destination file> <colors count>");
        return -1;
    }

    Image(argv[1], argv[2], std::stoi(argv[3]));

    return 0;
}
