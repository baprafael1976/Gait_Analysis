Análise de Forças na Marcha a partir de Arquivo .emt

Rotina em Python para ler arquivos de força no tempo (.emt) de duas plataformas de força em sequência, aplicar filtragem, calcular integrais de força, simetria, fases de suporte, velocidade média e estimativas de gasto mecânico, além de gerar gráficos e um relatório em CSV.

Estrutura dos dados de entrada

Arquivo texto com 8 colunas por linha, separadas por espaço:

Frame | Time | Fx1 | Fy1 | Fz1 | Fx2 | Fy2 | Fz2

Observações
	•	A versão original do código ignorava as duas primeiras linhas; a versão adaptada lê todas as linhas e descarta linhas onde Fz1 e Fz2 são ambos NaN.
	•	Amostragem assumida: 500 Hz.
	•	Unidades esperadas: força em N e tempo em s.

Principais funcionalidades
	•	Leitura e limpeza dos dados com descarte de linhas inválidas.
	•	Filtragem passa-baixa de Butterworth (10 Hz) em Fx, Fy, Fz de cada pé.
	•	Gráficos das componentes de força ao longo do tempo (X, Y, Z) e gráfico dedicado a Y.
	•	Índices de simetria por pico para Fx, Fy e Fz.
	•	Integrais de Fy por pé via regra do trapézio.
	•	Estimativas simples de deslocamento, trabalho e gasto mecânico a partir de integrais força-tempo.
	•	Estimativa de fases de suporte simples, duplo e contagem de passos a partir de NaN em Fy1/Fy2.
	•	Estimativas de tempo total, número de passos, cadência, tempo de passo, velocidade média e velocidade em km/h.
	•	Exporta um relatório acumulado em CSV.

Saídas
	•	Gráficos exibidos na tela com matplotlib.
	•	Relatório CSV em /content/report.csv contendo uma linha por arquivo analisado com as chaves abaixo.

Chaves do relatório e significado
	•	Total Time [s]
	•	Number of Steps
	•	Symmetry Index Fx, Symmetry Index Fy, Symmetry Index Fz
	•	Integral Fy (R) [N·s], Integral Fy (L) [N·s]
	•	Mechanical Energy Expenditure (R) [J], Mechanical Energy Expenditure (L) [J]
	•	Single Support Phase (R) [s], Single Support Phase (L) [s], Double Support Phase [s]
	•	Average Velocity [m/s], Average Velocity [km/h]
	•	Step Frequency [step/s]
	•	Step Time (Number of Steps Based) [s], Step Time (Support Phase Based) [s]

Requisitos
	•	Python 3.10+
	•	numpy, pandas, matplotlib, scipy

Instalação rápida

pip install numpy pandas matplotlib scipy

Como usar
	1.	Coloque o arquivo .emt no mesmo diretório do script.
	2.	Ajuste no final do script a variável file_path para o nome do seu arquivo, por exemplo:

file_path = "10-rapida-force.emt"
analyze_file(file_path)

	3.	Execute:

python seu_script.py

Notas sobre ambiente
	•	Se executar no Google Colab, o relatório será salvo em /content/report.csv. Em ambiente local, altere o caminho em export_report para um diretório de sua preferência.

Parâmetros e suposições
	•	Frequência de corte do filtro passa-baixa: 10 Hz.
	•	Frequência de amostragem: 500 Hz.
	•	Comprimento efetivo por plataforma: 0,60 m; número de plataformas: 4.
	•	Cálculo de velocidade média assume travessia de 4 plataformas em linha (2,40 m) e usa o tempo total do trecho analisado.

Limitações e pontos de atenção
	•	A rotina de gasto mecânico usa integrais de força para obter velocidade e deslocamento sem modelar massa ou cinemática real; interprete como estimativa grosseira e não como custo metabólico ou trabalho mecânico segmentado clássico.
	•	As fases de suporte são inferidas a partir de NaN em Fy1/Fy2; isso depende de como seu arquivo representa ausência de contato. Se sua aquisição não usa NaN, adapte a lógica de detecção de contato (por exemplo, limiar de Fz).
	•	Existem duas definições de analyze_file no código original; esta versão usa a segunda (adaptada). Recomenda-se manter apenas uma.
	•	A função butter_bandpass_filter está incompleta e não é utilizada. Pode ser removida ou finalizada.
	•	Pequenos erros de grafia nas chaves do relatório (velociry_ms/velociry_kmh). Mantidos para compatibilidade, mas recomenda-se corrigir.
	•	O caminho do CSV está fixo em /content; para uso fora do Colab, torne esse caminho configurável.

Checklist de reprodutibilidade
	•	Informe a versão do Python e das bibliotecas (pip freeze ou arquivo requirements.txt).
	•	Documente a origem do arquivo .emt, taxa de amostragem real e unidades.
	•	Fixe uma versão do código via Git tag e, se possível, gere um DOI (GitHub conectado ao Zenodo).
	•	Salve os gráficos em arquivos com nomes determinísticos se precisar anexá-los em relatório.

Estrutura sugerida do repositório

.
├── README.md
├── requirements.txt
├── scripts/
│   └── analyze_emt.py
├── data/
│   └── exemplo.emt
└── outputs/
    └── report.csv

Exemplo de citação

Vancouver
Baptista RR. Rotina de análise de forças em marcha [software]. Versão 1.0. 2025. Disponível em: DOI a ser gerado no Zenodo.

APA
Baptista, R. R. (2025). Rotina de análise de forças em marcha (Versão 1.0) [Software]. DOI a ser gerado no Zenodo.

ABNT
BAPTISTA, R. R. Rotina de análise de forças em marcha. Versão 1.0. 2025. DOI a ser gerado no Zenodo.

Licença

Defina uma licença explícita no repositório. Sugestões comuns para código científico: MIT ou Apache-2.0.

Roadmap curto
	•	Corrigir nomes das chaves velociry_ms/velociry_kmh.
	•	Tornar caminho do relatório configurável via argumento de linha de comando.
	•	Unificar analyze_file e remover funções não usadas.
	•	Tornar detectores de contato baseados em Fz com limiar configurável.
	•	Opcional: salvar figuras em outputs/ com nomes padronizados.

Contato

Rafael Reimann Baptista — PUCRS
Se utilizar este código em trabalhos acadêmicos, por favor cite a versão específica (de preferência com DOI).
