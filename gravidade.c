/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Algoritimo para simulação de efeitos gravitacionais em duas dimensões *
* Data: 03/04/2020                                                      *
* Autor: Vitor Henrique Andrade Helfensteller Straggiotti Silva         *
* FALTA CODIFICAR A SAIDA DE DADOS PARA UM ARQUIVO TXT                  *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

//estrutura que descreve o corpo puntiforme
struct corpo{
  float massa;
  float X, Y;     //coordenadas de posição
  float Vx, Vy;   //coordenadas de velocidade
  float Acx, Acy; //coordenadas de aceleração
  float Fx, Fy;   //coordenadas de Força
}; typedef struct corpo corpo;  //redefinindo escrita de estrutura
//estrutura que descreve as coordenadas de posição de um corpo
struct posicao{
  float X, Y;
}; typedef struct posicao posicao; //redefinindo escrita de estrutura

float    mod_cubo(corpo* objeto, int index_obj_1, int index_obj_2);
float**  matriz_colisao(corpo* objeto, int Num_Obj, float Dist_min);

int main(void){
  //Constantes
  const float G = 0.000000000066743; //Constante de gravitação universal
  //Variaveis de programa
  float delta_t = 0.001;                // (1ms) elemento de discretização do tempo
  unsigned int Itera_lim = 1000000;     // (1000s) numero de iterações
  unsigned int Janela_tmp = 10000;     // (10s) intervalo de tempo entre impressões do resultado
  float tempo = 0;                      // Armazena o tempo em que a simulação se encontra
  int Num_Obj = 2;                      // Numero de objetos a serem simulados
  float Dist_min = 0.25;                // Distancia minima entre as particulas "contato"

  corpo objeto[Num_Obj];                // Estruturas para armazenar variaveis do objeto
  posicao Posicao_anterior[Num_Obj];    // Estruturas para armazenar coordenadas anteriores de posição
  float** Matriz_col;                   // Matriz triangular superior que armazena estado binario de colisao

  //Inicializando variaveis dos corpos
  for(int i=0; i<Num_Obj; i++){
    objeto[i].massa = 0;
    objeto[i].X = 0;    objeto[i].Y = 0;
    objeto[i].Vx = 0;   objeto[i].Vy = 0;
    objeto[i].Acx = 0;  objeto[i].Acy = 0;
    objeto[i].Fx = 0;   objeto[i].Fy = 0;
  }

  //Adquirindo condições iniciais dos objetos
  for(int i=0; i<Num_Obj; i++){
    printf("Informe a massa do objeto %d : ", i);
    scanf("%f", &objeto[i].massa);
    printf("Informe a coordenada inicial de posicao X do objeto %d : ", i);
    scanf("%f", &objeto[i].X);
    printf("Informe a coordenada inicial de posicaoY do objeto %d : ", i);
    scanf("%f", &objeto[i].Y);
    printf("Informe a coordenada inicial de velocidade Vx do objeto %d : ", i);
    scanf("%f", &objeto[i].Vx);
    printf("Informe a coordenada inicial de velocidade Vy do objeto %d : ", i);
    scanf("%f", &objeto[i].Vy);
  }

  //Inicializando e declarando arquivo de saida
  FILE* Arquivo_saida;
  Arquivo_saida = fopen("simulation_out.txt", "w");
  if(Arquivo_saida == NULL){
    printf("ERRO!! Nao foi possivel criar arquivo.");
    exit(1);
  }//fim de if

  //Loop para calculo das variaveis cinematicas
  for(unsigned int Itera=0; Itera<Itera_lim; Itera++){
    //preenchimento da memoria de posição antes de iteração cinemática
    for(int i=0; i<Num_Obj; i++){
      Posicao_anterior[i].X = objeto[i].X;
      Posicao_anterior[i].Y = objeto[i].Y;
    }

    //Calculo das forças gavitacionais entre os objetos
    for(int i=0; i<Num_Obj; i++){
      for(int j=0; j<Num_Obj; j++){
        //para cada objeto i percorremos todos os objetos com j
        if(i != j){  //para evitar dividir por zero
          objeto[i].Fx = objeto[i].Fx + (G * objeto[i].massa * objeto[j].massa * (objeto[j].X - objeto[i].X))/(mod_cubo(objeto, i, j));
          objeto[i].Fy = objeto[i].Fy + (G * objeto[i].massa * objeto[j].massa * (objeto[j].Y - objeto[i].Y))/(mod_cubo(objeto, i, j));
        }//fim do if
      }//fim do for j
    }//fim do for i

    //Calculo das acelerações de cada objeto
    for(int i=0; i<Num_Obj; i++){
      objeto[i].Acx = objeto[i].Fx/objeto[i].massa;
      objeto[i].Acy = objeto[i].Fy/objeto[i].massa;
    }

    //Calculo de velocidades e posições
    for(int i=0; i<Num_Obj; i++){
      //Velocidades V=(X-Xo)/dt ==> X=Xo+Vdt
      objeto[i].Vx = objeto[i].Vx + (objeto[i].Acx * delta_t);
      objeto[i].Vy = objeto[i].Vy + (objeto[i].Acy * delta_t);
      //Posições a=(V-Vo)/dt ==> V=Vo+adt
      objeto[i].X = objeto[i].X + (objeto[i].Vx * delta_t);
      objeto[i].Y = objeto[i].Y + (objeto[i].Vy * delta_t);
    }

    //Testando condição de distancia mínima ("corpos em contato")
    Matriz_col = matriz_colisao(objeto, Num_Obj, Dist_min);
    for(int i=0; i<Num_Obj; i++){
      for(int j=0; j<Num_Obj; j++){
        //para cada objeto i varremos todos os objetos com j
        if((j>i) && (Matriz_col[i][j])){ //caso haja colisão, recuperar valores anteriores de posição
          objeto[i].X = Posicao_anterior[i].X;
          objeto[i].Y = Posicao_anterior[i].Y;
          objeto[j].X = Posicao_anterior[j].X;
          objeto[j].Y = Posicao_anterior[j].Y;
        }//fim do if
      }//fim do for j
    }//fim do for i
    //Liberação de memoria alocada
    for(int i=0; i<Num_Obj; i++){
      free(Matriz_col[i]);
    }
    free(Matriz_col);

    //Atualizando tempo da simulação
    tempo = Itera * delta_t;

    //exportando dados para arquivo texto
    //impressão do resultado
    if((Itera%Janela_tmp)==0){
      fprintf(Arquivo_saida, "%f:", tempo);
      for(int i=0; i<Num_Obj; i++){
        fprintf(Arquivo_saida, "%d:", i);
        fprintf(Arquivo_saida, "%f,", objeto[i].X);
        fprintf(Arquivo_saida, "%f:", objeto[i].Y);
        fprintf(Arquivo_saida, "%f,", objeto[i].Vx);
        fprintf(Arquivo_saida, "%f:", objeto[i].Vy);
        fprintf(Arquivo_saida, "%f,", objeto[i].Acx);
        fprintf(Arquivo_saida, "%f:", objeto[i].Acy);
        fprintf(Arquivo_saida, "%f,", objeto[i].Fx);
        fprintf(Arquivo_saida, "%f/", objeto[i].Fy);
      }//fim de for i
      fprintf(Arquivo_saida, "\n");
    }//fim de if

    //zerando forças antes de realizar nova simulação
    for(int i=0; i<Num_Obj; i++){
      objeto[i].Fx = 0;
      objeto[i].Fy = 0;
    }

  }//fim do for Itera

  fclose(Arquivo_saida); //fechando arquivo de saida
  return 0;
}
/****************************************************************************/
//retorna o calculo do cubo da distancia entre dois corpos
float mod_cubo(corpo* objeto, int index_obj_1, int index_obj_2){
  float Resultado;
  float Dif_X, Dif_Y;
  Dif_X = objeto[index_obj_2].X - objeto[index_obj_1].X;
  Dif_Y = objeto[index_obj_2].Y - objeto[index_obj_1].Y;
  Resultado = (sqrt((Dif_X * Dif_X)+(Dif_Y * Dif_Y)))*(sqrt((Dif_X * Dif_X)+(Dif_Y * Dif_Y)))*(sqrt((Dif_X * Dif_X)+(Dif_Y * Dif_Y)));
  return Resultado;
}
//====================================================================
//retorna uma matriz triangular superior de estados binarios (1 ==> colisão)
float** matriz_colisao(corpo* objeto, int Num_Obj, float Dist_min){
  float** colisao;
  float Delta_X, Delta_Y;
  //Alocando memoria "distancia[Num_Obj][Num_Obj]"
  colisao = (float**)malloc(Num_Obj * sizeof(float*));
  for(int k=0;k<Num_Obj; k++){
    colisao[k] = (float*)malloc(Num_Obj * sizeof(float));
  }
  //Calculando distancia
  for(int i=0; i<Num_Obj; i++){
    for(int j=0; j<Num_Obj; j++){
      //para cada objeto i varremos todos os objetos com j
      Delta_X = objeto[j].X - objeto[i].X;
      Delta_Y = objeto[j].Y - objeto[i].Y;
      if(j>i){ //eliminamos as distancias repetidas (i==>j é igual a j==>i)
        if((sqrt((Delta_X*Delta_X) + (Delta_Y*Delta_Y))) < Dist_min){
          colisao[i][j] = 1;
        }else{
          colisao[i][j] = 0;
        }
      }//fim do if
    }//fim do for j
  }//fim do for i
  return colisao;
}
//=====================================================================
