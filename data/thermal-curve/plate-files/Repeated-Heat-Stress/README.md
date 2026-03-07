# Raw data for testing the effects of repeated thermal stress on Pacific oyster (Crassostrea gigas) metabolic performance using resazurin assay

Experiment ran by Christina Zhang
Reciprocal treatment plates set up as follows:
Plate 1: received no heat followed by no heat
Plate 2: received no heat followed by heat
Plate 3: received heat followed by no heat
Plate 4: received heat followed by heat
Plate 3 & 4 are exposed to the designed temperature (33 °C, 35 °C, or 37 °C) for 24 hours, then removed. All four plates then undergo daily water changes for 7 days. After 7 days, Plates 2 & 4 are exposed to the same designed temperature for 24 hours. Resazurin data is collected immediately after the second round of heat exposure.

## Data Naming Convention

Three rounds of data were collected for each experimental temperature.

### Round 1
**Folder:** `XXC-healthy-CZ`  
First round of data under the **XX °C temperature treatment**. Oysters are in **healthy** conditions (define as no unusual mortalities under room temperature for 7 days period of time during preliminary testing). **CZ** stands for Christina Zhang initials.

File naming format:
`plateX_TX`

- **X** = plate number  
- **T** = time point  
  - `T0` = immediately after resazurin was added  
  - `T0.5` = 30 minutes after resazurin was added

---

### Round 2
**Folder:** `XXC-healthy-R2-CZ`  
Second round of data under the **XX °C temperature treatment**.

File naming format:
`plateX_TX_R2`

- Same structure as Round 1  
- `R2` indicates **second round of data collection**

---

### Round 3
**Folder:** `XXC-healthy-R3-CZ`  
Third round of data under the **XX °C temperature treatment**.

File naming format:
`plateX_TX_R3`

- Same structure as above  
- `R3` indicates **third round of data collection**
