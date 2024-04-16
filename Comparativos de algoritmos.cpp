#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <random>
using namespace std;
using namespace std::chrono;

// Função para gerar um vetor aleatório de tamanho n
vector<int> generateRandomVector(int n) {
    // Inicializando o gerador de números aleatórios
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 100);

    // Vetor para armazenar números aleatórios
    vector<int> v(n);

    // Preencher o vetor com números aleatórios
    for (int i = 0; i < n; ++i) {
        v[i] = dis(gen);
    }

    return v;
}

// Função para trocar dois elementos em um vetor
void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

// Função de partição para o Quick Sort
int partition(vector<int>& v, int low, int high) {
    int pivot = v[high];
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (v[j] < pivot) {
            i++;
            swap(&v[i], &v[j]);
        }
    }
    swap(&v[i + 1], &v[high]);
    return (i + 1);
}

// Função principal do Quick Sort
void quickSort(vector<int>& v, int low, int high) {
    if (low < high) {
        int pi = partition(v, low, high);
        quickSort(v, low, pi - 1);
        quickSort(v, pi + 1, high);
    }
}

// Função de merge do merge sort
void merge(vector<int> &v, int low, int high, int mid){
    vector<int> c(high - low + 1); // Vetor temporário para armazenar os elementos mesclados
    int i = low; // Índice inicial para a primeira metade do vetor
    int j = mid + 1; // Índice inicial para a segunda metade do vetor
    int k = 0; // Índice inicial para o vetor temporário

    // Mescla as duas metades do vetor em um vetor temporário
    while (i <= mid && j <= high) {
        if (v[i] <= v[j]) {
            c[k] = v[i];
            i++;
        } else {
            c[k] = v[j];
            j++;
        }
        k++;
    }

    // Copia os elementos restantes da primeira metade (se houver)
    while (i <= mid) {
        c[k] = v[i];
        i++;
        k++;
    }

    // Copia os elementos restantes da segunda metade (se houver)
    while (j <= high) {
        c[k] = v[j];
        j++;
        k++;
    }

    // Copia os elementos do vetor temporário de volta para o vetor original
    for (int m = 0; m < k; m++) {
        v[low + m] = c[m];
    }
}

// Função principal do Merge Sort
void mergeSort(vector<int> &v, int low, int high){
    if(low < high){
        int mid = (low + high) / 2;
        mergeSort(v, low, mid);
        mergeSort(v, mid + 1, high);
        merge(v, low, high, mid);
    }
}

// Implementação do Selection Sort
void selectionSort(vector<int> &v) {
    int n = v.size();
    int min;
    for(int i = 0; i < n-1; i++) {
        min = i;
        for(int j = i+1; j < n; j++) {
            if (v[j] < v[min]) {
                min = j;
            }
        }
        if(min != i) {
            swap(&v[min], &v[i]);
        }
    }
}

// Implementação do Insertion Sort
void insertionSort(vector<int> &v) {
    int n = v.size();
    int aux, j;
    for(int i = 1; i < n; i++) {
        aux = v[i];
        for(j = i; j > 0 && aux < v[j-1]; j--) {
            v[j] = v[j-1];
        }
        v[j] = aux;
    }
}

// Implementação do Bubble Sort
void bubbleSort(vector<int> &v) {
    int n = v.size();
    for(int i = 0; i < n-1; i++) {
        for(int j = 0; j < n-i-1; j++) {
            if (v[j] > v[j+1]) {
                swap(&v[j], &v[j+1]);
            }
        }
    }
}

// Implementação do Shell Sort
void shellSort(vector<int> &v) {
    int n = v.size();
    for(int t = n/2; t > 0; t /= 2) {
        for(int i = t; i < n; i += 1) {
            int temp = v[i];
            int j;
            for(j = i; j >= t && v[j-t] > temp; j -= t) {
                v[j] = v[j-t];
            }
            v[j] = temp;
        }
    }
}

int main() {
    
    int n = 10; // tamanho do vetor v

    // Gerar vetor aleatório
    vector<int> v = generateRandomVector(n);

    // Copiar o vetor aleatório para outros vetores para usar em diferentes algoritmos
    vector<int> v_selection = v;
    vector<int> v_insertion = v;
    vector<int> v_bubble = v;
    vector<int> v_shell = v;
    vector<int> v_merge = v;
    vector<int> v_quick = v;


 //auto start_selection = high_resolution_clock::now(); Marca o tempo de início da execução do algoritmo
 //auto end_selection = high_resolution_clock::now();  Marca o tempo de término do termino do algoritmo
    auto start_selection = high_resolution_clock::now();    
    selectionSort(v_selection);
    auto end_selection = high_resolution_clock::now(); 

    auto start_insertion = high_resolution_clock::now(); 
    insertionSort(v_insertion);
    auto end_insertion = high_resolution_clock::now(); 

    auto start_bubble = high_resolution_clock::now(); 
    bubbleSort(v_bubble);
    auto end_bubble = high_resolution_clock::now(); 

    auto start_shell = high_resolution_clock::now(); 
    shellSort(v_shell);
    auto end_shell = high_resolution_clock::now(); 

    auto start_merge = high_resolution_clock::now(); 
    mergeSort(v_merge, 0, n - 1);
    auto end_merge = high_resolution_clock::now(); 

    auto start_quick = high_resolution_clock::now(); 
    quickSort(v_quick, 0, n - 1);
    auto end_quick = high_resolution_clock::now(); 

    // Calcular os tempos totais de execução para cada algoritmo
    duration<double> time_taken_selection = duration_cast<duration<double>>(end_selection - start_selection);
    duration<double> time_taken_insertion = duration_cast<duration<double>>(end_insertion - start_insertion);
    duration<double> time_taken_bubble = duration_cast<duration<double>>(end_bubble - start_bubble);
    duration<double> time_taken_shell = duration_cast<duration<double>>(end_shell - start_shell);
    duration<double> time_taken_merge = duration_cast<duration<double>>(end_merge - start_merge);
    duration<double> time_taken_quick = duration_cast<duration<double>>(end_quick - start_quick);

    cout << "\nTempo total para ordenação com Selection Sort: " << fixed << setprecision(9) << time_taken_selection.count() << " seg" << endl;
    cout << "\nTempo total para ordenação com Insertion Sort: " << fixed << setprecision(9) << time_taken_insertion.count() << " seg" << endl;
    cout << "\nTempo total para ordenação com Bubble Sort: " << fixed << setprecision(9) << time_taken_bubble.count() << " seg" << endl;
    cout << "\nTempo total para ordenação com Shell Sort: " << fixed << setprecision(9) << time_taken_shell.count() << " seg" << endl;
    cout << "\nTempo total para ordenação com Merge Sort: " << fixed << setprecision(9) << time_taken_merge.count() << " seg" << endl;
    cout << "\nTempo total para ordenação com Quick Sort: " << fixed << setprecision(9) << time_taken_quick.count() << " seg" << endl;

    return 0;
}
