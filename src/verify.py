def verificar_matriz_soma_uniforme(caminho_arquivo, soma_esperada_linha, soma_esperada_coluna):
    """
    Lê uma matriz binária de um arquivo e verifica se a soma de 1s em cada linha e
    coluna corresponde a um valor único esperado para todas as linhas e colunas.

    Args:
        caminho_arquivo (str): O caminho para o arquivo .txt contendo a matriz.
        soma_esperada_linha (int): O valor de soma esperado para TODAS as linhas.
        soma_esperada_coluna (int): O valor de soma esperado para TODAS as colunas.
    """
    matriz = []
    # Lê o arquivo e constrói a matriz
    try:
        with open(caminho_arquivo, 'r') as f:
            for linha in f:
                # Converte cada elemento da linha para int e adiciona à matriz
                matriz.append([int(x) for x in linha.strip().split()])
    except FileNotFoundError:
        print(f"Erro: O arquivo '{caminho_arquivo}' não foi encontrado.")
        return
    except ValueError:
        print("Erro: O arquivo contém caracteres que não são números inteiros.")
        return

    if not matriz:
        print("A matriz está vazia.")
        return

    num_linhas = len(matriz)
    # Assume que todas as linhas têm o mesmo número de colunas
    num_colunas = len(matriz[0])

    print("--- Verificando Linhas ---")
    check = True
    # Itera sobre cada linha para verificar a soma
    for i, linha in enumerate(matriz):
        soma_linha_atual = sum(linha)
        if soma_linha_atual != soma_esperada_linha:
            print(f"Linha {i}: Esperado: {soma_esperada_linha}, Encontrado: {soma_linha_atual}")
            check = False
    if check == True:
        print("Válido!")
            
    check = True
    print("\n--- Verificando Colunas ---")
    # Itera sobre cada coluna para verificar a soma
    for j in range(num_colunas):
        soma_coluna_atual = sum(matriz[i][j] for i in range(num_linhas))
        if soma_coluna_atual != soma_esperada_coluna:
            print(f"Coluna {j}: Esperado: {soma_esperada_coluna}, Encontrado: {soma_coluna_atual}")
            check = False
    if check == True:
        print("Válido!")


caminho_do_arquivo = 'cff_matrix.txt'
soma_esperada_para_cada_linha = 16
soma_esperada_para_cada_coluna = 4

print(f"Verificando a matriz com soma esperada de '{soma_esperada_para_cada_linha}' para as linhas e '{soma_esperada_para_cada_coluna}' para as colunas.\n")

verificar_matriz_soma_uniforme(
    caminho_do_arquivo,
    soma_esperada_para_cada_linha,
    soma_esperada_para_cada_coluna
)