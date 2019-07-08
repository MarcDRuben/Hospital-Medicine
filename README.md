Tools to study round-the-clock hospital medicine
================
Marc Ruben
7/8/2019

In a study <https://www.biorxiv.org/content/10.1101/617944v1> we
analyzed the daily distribution of ~120K doses of 12 separate drugs in
~1.5K inpatients at a major children’s hospital in the U.S. Treatment
orders and administration (first-doses) were strongly
time-of-day-dependent. These 24 h rhythms were consistent across drugs,
diagnoses, and hospital units.

![image caption Source](images/GitRepo_AllDrugWheels.png)

**This repository contains the code used for each part of this study.**

## Part 1. Extract EMR data from Epic

Extract patient EMR data from Epic. Requires SQL. Code may need to be
adapted to institutional differences in Epic.

## Part 2. Evaluate 24 h patterns in TREATMENT

Evaluate the daily distribution of treatment orders and first-doses
administered. Requires R.

## Part 3. Evaluate 24 h patterns in RESPONSE

Characterize the clinical response to hydralazine, an acutely
administered antihypertensive, as a function of time of day. Requires R.
