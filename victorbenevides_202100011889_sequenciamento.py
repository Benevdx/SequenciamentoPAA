import time
from dataclasses import dataclass
from typing import List
import multiprocessing


@dataclass
class DNA:
    tamanho_minimo_substring: int
    sequencia_dna: str


@dataclass
class Doenca:
    codigo: str
    genes: List[str]
    probabilidade_doenca: int = 0


def encontrar_ocorrencias_gene(sequencia_dna: str, gene: str, tamanho_minimo_substring: int) -> bool:
    """
    Correspondência de padrões super otimizada - apenas verifica se o gene aparece vezes suficientes
    para contribuir para a probabilidade da doença. Retorna antecipadamente quando o limite é atingido.
    """
    if not gene or not sequencia_dna:
        return False
        
    # Só precisa encontrar se o gene aparece vezes suficientes para contar
    contagem = 0
    inicio = 0
    
    # Busca de substring simples e rápida que retorna assim que o limite é atingido
    while inicio <= len(sequencia_dna) - len(gene):
        posicao = sequencia_dna.find(gene, inicio)
        if posicao == -1:
            break
            
        contagem += len(gene)
        if contagem >= tamanho_minimo_substring:
            return True
            
        inicio = posicao + 1
            
    return False


def calcular_probabilidade_doenca(sequencia_dna: str, genes: List[str], tamanho_minimo_substring: int):
    """Cálculo de probabilidade de doença ultra-otimizado"""
    if not genes:
        return 0
    
    # Usa uma abordagem de bitmap para melhor desempenho em listas grandes
    genes_correspondentes = 0
    
    for gene in genes:
        if encontrar_ocorrencias_gene(sequencia_dna, gene, tamanho_minimo_substring):
            genes_correspondentes += 1
            
    probabilidade_doenca = round((genes_correspondentes / len(genes)) * 100)
    return min(probabilidade_doenca, 100)


def processar_doenca(argumentos):
    """Processa uma única doença e retorna o resultado"""
    linha_doenca, sequencia_dna, tamanho_minimo_substring = argumentos
    
    partes = linha_doenca.strip().split()
    codigo = partes[0]
    genes = partes[2:]  # Genes começam no terceiro elemento
    
    # Calcula a probabilidade da doença
    probabilidade_doenca = calcular_probabilidade_doenca(sequencia_dna, genes, tamanho_minimo_substring)
    
    return Doenca(codigo, genes, probabilidade_doenca)


def ler_arquivo(nome_arquivo_entrada: str):
    """Lê o arquivo de entrada com mínimo uso de memória"""
    with open(nome_arquivo_entrada, "r") as arquivo:
        tamanho_minimo_substring = int(arquivo.readline().strip())
        sequencia_dna = arquivo.readline().strip()
        quantidade_doencas = int(arquivo.readline().strip())
        linhas_doencas = [arquivo.readline().strip() for _ in range(quantidade_doencas)]
    
    return tamanho_minimo_substring, sequencia_dna, linhas_doencas


def ordenar_doencas(doencas: List[Doenca]):
    """
    Ordena doenças apenas por probabilidade (decrescente).
    Para mesmas probabilidades, mantém a ordem original (estabilidade).
    
    Implementação utilizando algoritmo de ordenação estável.
    """
    # Usando o algoritmo de ordenação estável do Python (Timsort)
    # Que mantém a ordem original para elementos com a mesma chave
    return sorted(doencas, key=lambda doenca: -doenca.probabilidade_doenca)


def escrever_arquivo(nome_arquivo_saida: str, doencas: List[Doenca]):
    """Escreve resultados em uma única operação"""
    with open(nome_arquivo_saida, "w") as arquivo:
        arquivo.writelines(f"{doenca.codigo}->{doenca.probabilidade_doenca}%\n" for doenca in doencas)


def main():
    nome_arquivo_entrada = "entrada.txt"
    nome_arquivo_saida = "saida.txt"

    # Medição de tempo
    tempo_inicio = time.time()

    # Obtém o número de núcleos de CPU disponíveis
    numero_nucleos = multiprocessing.cpu_count()
    
    # Lê dados do arquivo
    tamanho_minimo_substring, sequencia_dna, linhas_doencas = ler_arquivo(nome_arquivo_entrada)
    
    # Prepara argumentos para processamento paralelo
    argumentos_processamento = [(linha, sequencia_dna, tamanho_minimo_substring) for linha in linhas_doencas]
    
    # Usa pool de processos para paralelismo verdadeiro (evitando GIL)
    with multiprocessing.Pool(processes=numero_nucleos) as pool:
        # Mapeia o processamento de doenças para o pool
        doencas = pool.map(processar_doenca, argumentos_processamento)
    
    # Ordena doenças
    doencas_ordenadas = ordenar_doencas(doencas)

    # Escreve resultados
    escrever_arquivo(nome_arquivo_saida, doencas_ordenadas)

    # Tempo de execução
    tempo_fim = time.time()
    print(f"Tempo de execução: {tempo_fim - tempo_inicio:.6f} segundos")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()