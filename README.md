# Проект по майнору Биоинформатика, ВШЭ 2021

Целью работы над проектом является поиск и изучение участков генома, где
определенная гистоновая метка присутствует в местах образования одной из вторичных структур ДНК (ZDNA или квадруплекс, G4).


* **Организм:** human
* **Сборка генома:** hg19
* **Структура:** G4_seq_Li_KPDS
* **Метка:** H3K9me3
* **Тип:** A549

## Команды, используемые для выполнения работы

```python
# Скачивание файлов
wget https://www.encodeproject.org/files/ENCFF164FDB/@@download/ENCFF164FDB.bed.gz
wget https://www.encodeproject.org/files/ENCFF494QKI/@@download/ENCFF494QKI.bed.gz

# Расспаковка файлов с удалением ненужных столбцов
zcat ENCFF164FDB.bed.gz  |  cut -f1-5 > H3K9me3_A549.ENCFF164FDB.hg38.bed
zcat ENCFF494QKI.bed.gz  |  cut -f1-5 > H3K9me3_A549.ENCFF494QKI.hg38.bed

# Конвертация hg38 => hg19
liftOver   H3K9me3_A549.ENCFF164FDB.hg38.bed   hg38ToHg19.over.chain.gz   H3K9me3_A549.ENCFF164FDB.hg19.bed   H3K9me3_A549.ENCFF164FDB.unmapped.bed
liftOver   H3K9me3_A549.ENCFF494QKI.hg38.bed   hg38ToHg19.over.chain.gz   H3K9me3_A549.ENCFF494QKI.hg19.bed   H3K9me3_A549.ENCFF494QKI.unmapped.bed

# Скачивание скрипта для построения гистограм
wget https://github.com/vanya-antonov/hse21_H3K4me3_ZDNA_human/raw/main/src/lib.R
wget https://github.com/vanya-antonov/hse21_H3K4me3_ZDNA_human/raw/main/src/len_hist.R

# Настройка пользователя github
git config --global user.name "ARATOMI"
git config --global user.email "*hiden*"
git config --global core.autocrlf input
git config --global color.ui auto

# Копирование репозитория
git clone https://github.com/ARATOMI/hse21_H3K9me3_G4_human.git

# Скачивание скрипта для фильтрации пиков
wget https://github.com/vanya-antonov/hse21_H3K4me3_ZDNA_human/raw/main/src/filter_peaks.R

# Объединение двух наборов отфильтрованных пиков
cat  *.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K9me3_A549.merge.hg19.bed

# Скачивание Homo_Li_KPDS
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3003540&format=file&file=GSM3003540_Homo_all_w15_th-1_minus%2Ehits%2Emax%2EPDS%2Ew50%2E35%2Ebed%2Egz
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3003540&format=file&file=GSM3003540_Human_all_w15_th-1_minus%2Ecov%2EPDS%2EbedGraph%2Egz

# Распаковка полученных файлов
gunzip GSM3003540_Homo_all_w15_th-1_minus.hits.max.PDS.w50.35.bed.gz
gunzip GSM3003540_Homo_all_w15_th-1_plus.hits.max.PDS.w50.35.bed.gz

# Объединение в 1 файл
cat  *.hits.max.PDS.w50.35.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  G4_Li_KPDS.bed

# Команды для Genome Browser
track visibility=dense name="ENCFF164FDB"  description="3K9me3_A549.ENCFF164FDB.hg19.filtered.bed"
https://raw.githubusercontent.com/ARATOMI/hse21_H3K9me3_G4_human/main/data/H3K9me3_A549.ENCFF164FDB.hg19.filtered.bed

track visibility=dense name="ENCFF494QKI"  description="H3K9me3_A549.ENCFF494QKI.hg19.filtered.bed"
https://raw.githubusercontent.com/ARATOMI/hse21_H3K9me3_G4_human/main/data/H3K9me3_A549.ENCFF494QKI.hg19.filtered.bed

track visibility=dense name="ChIP_merge"  color=50,50,200   description="H3K9me3_A549.merge.hg19.bed"
https://raw.githubusercontent.com/ARATOMI/hse21_H3K9me3_G4_human/main/data/H3K9me3_A549.merge.hg19.bed

track visibility=dense name="G4_Li_KPDS"  color=0,200,0   description="G4_Li_KPDS.bed"
https://raw.githubusercontent.com/ARATOMI/hse21_H3K9me3_G4_human/main/data/G4_Li_KPDS.bed

# Нахождение пересечение гистоновой метки и структуры ДНК
bedtools intersect  -a G4_Li_KPDS.bed   -b  H3K9me3_A549.merge.hg19.bed  >  H3K9me3_A549.intersect_with_G4_Li_KPDS.bed

# Команды для Genome Browser
track visibility=dense name="intersect_with_G4_Li_KPDS"  color=200,0,0  description="H3K9me3_A549.intersect_with_G4_Li_KPDS.bed"
https://raw.githubusercontent.com/ARATOMI/hse21_H3K9me3_G4_human/main/data/H3K9me3_A549.intersect_with_G4_Li_KPDS.bed

# Скачивание скрипта ChIPpeakAnno.R
wget https://raw.githubusercontent.com/vanya-antonov/hse21_H3K4me3_ZDNA_human/main/src/ChIPpeakAnno.R


```

## Результаты работы

![len_hist.H3K9me3_A549.ENCFF494QKI.hg38.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/len_hist.H3K9me3_A549.ENCFF494QKI.hg38.png)

![len_hist.H3K9me3_A549.ENCFF164FDB.hg38.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/len_hist.H3K9me3_A549.ENCFF164FDB.hg38.png)

![len_hist.H3K9me3_A549.ENCFF494QKI.hg19.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/len_hist.H3K9me3_A549.ENCFF494QKI.hg19.png)

![len_hist.H3K9me3_A549.ENCFF164FDB.hg19.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/len_hist.H3K9me3_A549.ENCFF164FDB.hg19.png)

![len_hist.G4_Li_KPDS.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/len_hist.G4_Li_KPDS.png)

![len_hist.H3K9me3_A549.intersect_with_G4_Li_KPDS.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/len_hist.H3K9me3_A549.intersect_with_G4_Li_KPDS.png)

![chip_seeker.H3K9me3_A549.ENCFF494QKI.hg19.filtered.plotAnnoPie.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/chip_seeker.H3K9me3_A549.ENCFF494QKI.hg19.filtered.plotAnnoPie.png)

![chip_seeker.H3K9me3_A549.ENCFF164FDB.hg19.filtered.plotAnnoPie.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/chip_seeker.H3K9me3_A549.ENCFF164FDB.hg19.filtered.plotAnnoPie.png)

![chip_seeker.G4_Li_KPDS.plotAnnoPie.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/chip_seeker.G4_Li_KPDS.plotAnnoPie.png)

![chip_seeker.H3K9me3_A549.intersect_with_G4_Li_KPDS.plotAnnoPie.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/chip_seeker.H3K9me3_A549.intersect_with_G4_Li_KPDS.plotAnnoPie.png)

![gb.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/screenshots/gb.png)

![gb2.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/screenshots/gb2.png)

![go.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/screenshots/go.png)

![peaks.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/screenshots/peaks.png)

![wc.png](https://github.com/ARATOMI/hse21_H3K9me3_G4_human/raw/main/images/screenshots/wc.png)
