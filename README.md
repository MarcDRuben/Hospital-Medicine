Tools to study round-the-clock hospital medicine
================
Marc Ruben
7/8/2019

**This repository describes (with code) how we leveraged electronic
medical record (EMR) data to study treatment and response as a function
of time-of-day**

## Background

We analyzed the daily distribution of ~120K doses of 12 separate drugs
in ~1.5K inpatients at a major children’s hospital in the U.S. Treatment
orders and first-doses administered were strongly time-of-day-dependent
(figure below). These 24 h rhythms were consistent across drugs,
diagnoses, and hospital units.
<https://www.biorxiv.org/content/10.1101/617944v1>

![image caption Source](images/GitRepo_AllDrugWheels.png)

This repository contains descriptions (and code) to reproduce these
analyses.

## Part 1. Extract EMR data from Epic

Extract patient EMR data from Epic. Requires SQL. Code will need to be
adapted to institutional differences in how Epic data is maintained.

## Part 2. Evaluate 24 h patterns in TREATMENT

Characterize the daily distribution of treatment orders and first-doses
administered. Requires R.

## Part 3. Evaluate 24 h patterns in RESPONSE

Characterize clinical responses as a function of time of day. Requires
R.
