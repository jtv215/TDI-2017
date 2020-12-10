
#include <stdio.h>
#include <C_General.hpp>
#include <C_Trace.hpp>
#include <C_File.hpp>
#include <C_Arguments.hpp>
#include <C_Matrix.hpp>
#include <C_Image.hpp>
#include <iostream>
#include <cmath>
#include <conio.h>
#include <string>
#include <limits>

/*****VARIABLES******/
//Variables de posicionamiento y tamaño en la matriz de imagen a tratar
int firstRow, lastRow, firstCol, lastCol, nRows, nCols;
//Imagenes
C_Image imagenOriginal, imagenGauss, imagenLaplaciana,  imagenLaplaceInvertido, imagenLaplacianoAbsoluto, imagenLaplacianoInvertida, imagenPasoPorCero, imagenFinal;

//Inicializa las variables
void init(int mainMaxUmbral) {
	//Limites de la matriz imagen original
	firstRow = imagenOriginal.FirstRow();
	lastRow = imagenOriginal.LastRow();
	firstCol = imagenOriginal.FirstCol();
	lastCol = imagenOriginal.LastCol();
	//Tamaño de la matriz imagen original
	nRows = imagenOriginal.RowN();
	nCols = imagenOriginal.ColN();
}

//mascara Gaussiana de 3x3
C_Matrix mascaraGauss3x3() {
	C_Matrix filtro = C_Matrix(0, 2, 0, 2, 0);
	//Rellenamos la matriz de la siguiente forma:
	//___________
	//|_1_|_2_|_1_|
	//|_2_|_4_|_2_|
	//|_1_|_2_|_1_|
	filtro(0, 0) = 1; filtro(0, 1) = 2; filtro(0, 2) = 1;
	filtro(1, 0) = 2; filtro(1, 1) =4; filtro(1, 2) = 2;
	filtro(2, 0) = 1; filtro(2, 1) = 2; filtro(2, 2) = 1;
	return filtro;
}

//Mascara Gaussiana de 5x5
C_Matrix mascaraGauss5x5() {
	//________________________
	//|__1_|_4_|_6_|_4_|__1_|
	//|_4_|_16_|_24_|_16_|_4_|
	//|_6_|_24_|_36_|_24_|_6_|
	//|_4_|_16_|_24_|_16_|_4_|
	//|__1_|_4_|_6_|_4_|__1_|
	C_Matrix filtro = C_Matrix(0, 4, 0, 4, 0);
	filtro(0, 0) = 1; filtro(0, 1) = 4; filtro(0, 2) = 6; filtro(0, 3) = 4; filtro(0, 4) = 1;
	filtro(1, 0) = 4; filtro(1, 1) = 16; filtro(1, 2) = 24; filtro(1, 3) = 16; filtro(1, 4) = 4;
	filtro(2, 0) = 6; filtro(2, 1) = 24; filtro(2, 2) = 36; filtro(2, 3) = 24; filtro(2, 4) = 6;
	filtro(3, 0) = 4; filtro(3, 1) = 16; filtro(3, 2) = 24; filtro(3, 3) = 16; filtro(3, 4) = 4;
	filtro(4, 0) = 1; filtro(4, 1) = 4; filtro(4, 2) = 6; filtro(4, 3) = 4; filtro(4, 4) = 1;
	return filtro;
}




//Genera el filtro Laplaciano de 3x3
C_Matrix mascaraLaplace3x3() {
	C_Matrix filtro = C_Matrix(0, 2, 0, 2, 0);
	//____________
	//|_0_|_-1_|_0_|
	//|_-1_|_4_|_-1_|
	//|_0_|_-1_|_0_|
	filtro(0, 0) = 0; filtro(0, 1) = -1 ; filtro(0, 2) = 0;
	filtro(1, 0) = -1; filtro(1, 1) = 4 ; filtro(1, 2) = -1;
	filtro(2, 0) = 0; filtro(2, 1) = -1 ; filtro(2, 2) = 0;
	return filtro;
}
//Genera el filtro Laplaciano de 5x5
C_Matrix mascaraLaplace5x5() {
	C_Matrix filtro = C_Matrix(0, 4, 0, 4, 0);
	//________________________
	//|__0_|__0_|_-1_|__0_|__0_|
	//|_0_|_-1_|_-2_|_-1_|__0_|
	//|_-1_|_-2_|_16_|_-2_|_-1_|
	//|__0_|_-1_|_-2_|_-1_|_0_|
	//|__0_|__0_|_-1_|_0_|__0_|
	filtro(0, 0) = 0; filtro(0, 1) = 0; filtro(0, 2) = -1; filtro(0, 3) = 0; filtro(0, 4) = 0;
	filtro(1, 0) = 0; filtro(1, 1) = -1; filtro(1, 2) = -2; filtro(1, 3) = -1; filtro(1, 4) = 0;
	filtro(2, 0) = -1; filtro(2, 1) = -2; filtro(2, 2) = 16; filtro(2, 3) = -2; filtro(2, 4) = -1;
	filtro(3, 0) = 0; filtro(3, 1) = -1; filtro(3, 2) = -2; filtro(3, 3) = -1; filtro(3, 4) = 0;
	filtro(4, 0) = 0; filtro(4, 1) = 0; filtro(4, 2) = -1; filtro(4, 3) = 0; filtro(4, 4) = 0;
	return filtro;
}



C_Image convolucionValoresAbsolutos(C_Image imagen, C_Matrix filtro) {
	/*Generamos una matriz de salida copia de la imagen Gaussiana
	y establecemos todos sus valores a 0 (negro)*/
	C_Matrix out = imagen;
	out.SetValue(0);
	//Variables tamaño del filtro
	int fSizeX = filtro.ColN();
	int fSizeY = filtro.RowN();
	//Variables centro del filtro
	int fCenterX = fSizeX / 2;
	int fCenterY = fSizeY / 2;
	//Variable sumatoria
	int sum;
	//Variables auxiliares
	int rowIndex, colIndex;
	//Recorremos la matriz imagen
	for (int i = firstRow; i<nRows; i++) {
		for (int j = firstCol; j<nCols; j++) {
			//Inicializamos la sumatoria a 0
			sum = 0;
			//Recorremos la matriz del filtro
			for (int m = 0; m<fSizeY; m++) {
				for (int n = 0; n<fSizeX; n++) {
					rowIndex = i + m - fCenterY;
					colIndex = j + n - fCenterX;
					if (rowIndex >= firstRow && rowIndex<nRows && colIndex >= firstCol && colIndex<nCols) {
						//pixel de la imagen
						int a = imagen(rowIndex, colIndex);
						//pixel del flitro
						int b = filtro(m, n);
						//se multiplican y se unen a la suma
						sum += a*b;
						//En caso de que estemos cerca del borde, hacemos que
						//esos pixeles sean iguales que 0 (negro) y comprobamos
						//el siguiente pixel, es decir, salimos de la
						//comprobación del pixel actual
					}
					else {
						//Hacemos que la suma sea igual a un pixel 0 (negro)
						sum = 0;
						//Ponemos los indices m y n al tamaño del filtro
						//para salir del bucle
						m = fSizeY;
						n = fSizeX;
					}
				}
			}
			//Se introduce para ese pixel el resultado del computo de la convolución
			//en la matriz de salida
			out(i, j) = abs(sum);
		}
	}
	//Si la imagen excede el rango de 255, hacemos que la imagen se ajuste sus valores
	// para que se ajuste al rango [0,255]
	double max = out.Max();
	if (max>255.0) {
		out.Stretch(0, 255);
	}
	return out;
}


C_Matrix convolucion(C_Image imagen, C_Matrix filtro) {
	/*Generamos una matriz de salida copia de la imagen Gaussiana
	y establecemos todos sus valores a 0 (negro)*/
	C_Matrix out = imagen;
	out.SetValue(0);
	//Variables tamaño del filtro
	int fSizeX = filtro.ColN(); 
	int fSizeY = filtro.RowN();
	//Variables centro del filtro
	int fCenterX = fSizeX / 2;
	int fCenterY = fSizeY / 2;
	//Variable sumatoria
	int sum;
	//Variables auxiliares
	int rowIndex, colIndex;
	//Recorremos la matriz imagen
	for (int i = firstRow; i<nRows; i++) {
		for (int j = firstCol; j<nCols; j++) {
			//Inicializamos la sumatoria a 0
			sum = 0;
			//Recorremos la matriz del filtro
			for (int m = 0; m<fSizeY; m++) {
				for (int n = 0; n<fSizeX; n++) {
					rowIndex = i + m - fCenterY;
					colIndex = j + n - fCenterX;
					if (rowIndex >= firstRow && rowIndex<nRows && colIndex >= firstCol && colIndex<nCols) {
						//pixel de la imagen
						int a = imagen(rowIndex, colIndex);
						//pixel del flitro
						int b = filtro(m, n);
						//se multiplican y se unen a la suma
						sum += a*b;
						//En caso de que estemos cerca del borde, hacemos que
						//esos pixeles sean iguales que 0 (negro) y comprobamos
						//el siguiente pixel, es decir, salimos de la
						//comprobación del pixel actual
					}
					else {
						//Hacemos que la suma sea igual a un pixel 0 (negro)
						sum = 0;
						//Ponemos los indices m y n al tamaño del filtro
						//para salir del bucle
						m = fSizeY;
						n = fSizeX;
					}
				}
			}
			//Se introduce para ese pixel el resultado del computo de la convolución
			//en la matriz de salida
			out(i, j) = (sum);
		}
	}
	//Si la imagen excede el rango de 255, hacemos que la imagen se ajuste sus valores
	// para que se ajuste al rango [0,255]
	double max = out.Max();
	if (max>255.0) {
		out.Stretch(0, 255);
	}
	return out;
}







C_Image LaplacianoInvertido(C_Image imagenLaplacianoAbsoluto, int umbralMax) {
	/*Creamos una imagen que es la copia*/
	C_Image  imagenAux= imagenLaplacianoAbsoluto;
	/*Recorremos la imagen original*/
	for (int x = firstRow; x < lastRow; x++) {
		for (int y = firstCol; y < lastCol; y++) {
			if (imagenLaplacianoAbsoluto(x, y) < umbralMax)
				imagenAux(x, y) = 255;
			else
				imagenAux(x, y) = 0;
		}
	}
	/*Si nuestra imagen excede del rango 255, hacemos que la imagen se ajuste*/
	double max = imagenAux.Max();
	if (max>255.0) {
		imagenAux.Stretch(0, 255);
	}
	return imagenAux;
}


C_Image pasoPorCero(C_Image imagenLaplaciana, int umbralMax) {
	//____________
	//|p1_|p2_|p3_|
	//|p4_|_0_|p5_|
	//|p6_|p7_|p8_|
	double p0, p1, p2, p3, p4, p5, p6, p7, p8;

	//Generamos una imagen de tamaño la imagen del gradiente y establecemos sus valores a 0
	C_Image aux = imagenLaplaciana;
	aux.SetValue(250);
	//Recorremos la imagen de la place
	for (int i = firstRow + 1; i<lastRow; i++) {
		for (int j = firstCol + 1; j<lastCol; j++) {
			p0 = imagenLaplaciana(i, j);
			//Extraemos en las variables auxiliares los pixeles del
			//laplaciano de Derecha e Izquierda
			// ___________
			//|___|___|___|
			//|p4_|p0_|p5_|
			//|___|___|___|
			p4 = imagenLaplaciana(i, j - 1);
			p5 = imagenLaplaciana(i, j + 1);
			//Extraemos en las variables auxiliares los pixeles del
			//laplaciano de Arriba y abajo
			// ___________
			//|___|p2_|___|
			//|___|p0_|___|
			//|___|p7_|___|
			p2 = imagenLaplaciana(i - 1, j);
			p7 = imagenLaplaciana(i + 1, j);
			//Extraemos en las variables auxiliares los pixeles del
			//laplaciano de la diagonal1
			// ___________
			//|p1_|___|___|
			//|___|p0_|___|
			//|___|___|p8_|
			p1 = imagenLaplaciana(i - 1, j - 1);
			p8 = imagenLaplaciana(i + 1, j + 1);
			//Extraemos en las variables auxiliares los pixeles del
			//laplaciano de la diagonal2
			// ___________
			//|___|___|p3_|
			//|___|p0_|___|
			//|p6_|___|___|
			p3 = imagenLaplaciana(i - 1, j + 1);
			p6 = imagenLaplaciana(i + 1, j - 1);
			/*Los píxeles del borde son aquellos tal que el Laplaciano de dos de sus
			vecinos en posiciones opuestas tienen distinto signo
			Normalmente se considera un valor umbral para el valor absoluto de la
			diferencia numérica entre posiciones opuestas para considerar que un píxel
			es de paso por cero.*/
			if (abs(p4) - abs(p5) > umbralMax) {
				aux(i, j) = 0;//negro
			}
			else if (abs(p2) - abs(p7) > umbralMax) {
				aux(i, j) = 0;//negro
			}
			else if (abs(p1) - abs(p8) > umbralMax) {
				aux(i, j) = 0;//negro
			}
			else if (abs(p3) - abs(p6) > umbralMax) {
				aux(i, j) = 0;//negro
			}
			else {
				aux(i, j) == 255;//blanco
			}


		}
	}
	/*Si nuestra imagen excede del rango 255, hacemos que la imagen se
	ajuste*/
	double max = aux.Max();
	if (max>255.0) {
		aux.Stretch(0, 255);
	}
	return aux;
}

C_Image realzarBorde(C_Image imagen, C_Image imagenPasoCero, int umbral) {
	//Creamos una copia de la imagen original
	C_Image aux = imagen;
	//Modificamos la paleta de la imagen haciendo que el valor 0 corresponda a un color verde	
	aux.palette(0, C_GREEN) = 255;	
	//Recorremos la imagen de los bordes el paso por cero y si encontramos un pixel con valor 255, hacemos
	//que esa misma posicion valga verde en la imagen final
	for (int i = imagen.FirstRow(); i< imagen.LastRow(); i++) {
		for (int j = imagen.FirstCol(); j<imagen.LastCol(); j++) {	
			if (imagenPasoCero(i, j) < umbral) {
				aux(i, j) = 0;
			}
		}
	}
	return aux;
}

int main(int argc, char **argv)
{
	int umbralMax;
	char nombreImagen[50];
	cout << "Escribe el nombre de la imagen y su formato .bmp:" << endl;	
	cin.getline(nombreImagen, 50);
	imagenOriginal.ReadBMP(nombreImagen);

	if (imagenOriginal.Fail()) {
		cout << "La imagen no se encuentra. El programa se cerrara\n";
		system("pause");
		exit;

	}
	else {
		do {
			cout << "Indique el Umbral Maximo hasta de 1 al 100:\n";
			cin >> umbralMax;
			if (umbralMax <= 0 || umbralMax > 100) {
				cout << "Incorrecto, los valores deben de seres entre 1 a 100\n";
			}
		} while (umbralMax < 0 || umbralMax > 100);

		
		init(umbralMax);


		////////0.Imagen Original
		imagenOriginal.WriteBMP("0._ImagenOriginal.bmp");

		////////1.Gaus
		cout << "\n--------------------------------------\n";
		C_Print("Realizando tecnicas de detencion de bordes laplaciano");
		cout << "\n1.Tecnica de filtrado Gaussiano.\n";
		imagenGauss = convolucion(imagenOriginal, mascaraGauss5x5());
		imagenGauss.WriteBMP("1.imagen_Gaussiana.bmp");
				
		///////2.Laplaciano				
		cout << "\n2.Tecnica de filtrado Laplaciano.\n";
		imagenLaplaciana = convolucion(imagenGauss, mascaraLaplace3x3());
		imagenLaplaciana.WriteBMP("2.0_imagen_Laplaciana.bmp");

		cout << "\n 2.1Tecnica de Laplaciano con valores absolutos.\n";
		imagenLaplacianoAbsoluto = convolucionValoresAbsolutos(imagenGauss, mascaraLaplace3x3());
		imagenLaplacianoAbsoluto.WriteBMP("2.1_imagen_LaplacianoAbsolutos.bmp");
		
		cout << "\n 2.2Tecnica de Laplaciano Invertido\n";
		imagenLaplacianoInvertida = LaplacianoInvertido(imagenLaplacianoAbsoluto, umbralMax);
		imagenLaplacianoInvertida.WriteBMP("2.2_imagen_LaplacianoInvertida.bmp");
		
		///////3.Paso por cero
		cout << "\n3.Tecnica Paso por cero.\n";
		imagenPasoPorCero = pasoPorCero(imagenLaplaciana, umbralMax);
		imagenPasoPorCero.WriteBMP("3.Imagen_PasoPorCero.bmp");
			
		///////Resaltar los bordes
		cout << "\n4.Resalto de bordes.\n";
		imagenFinal = realzarBorde(imagenOriginal, imagenPasoPorCero, umbralMax);
		imagenFinal.WriteBMP("4.imagen_realzeDeBordes_Final.bmp");
		cout << "\n--------------------------------------\n";
				
		cout << "\nEl proceso se ha realizado con exito\n";
		cout << "\nAutor: Jefferson Max T.\n";


		system("pause");
		exit;
	}
	
}