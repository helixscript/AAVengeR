-- MariaDB dump 10.19  Distrib 10.6.12-MariaDB, for debian-linux-gnu (x86_64)
--
-- Host: 174.129.238.44    Database: AAVengeR
-- ------------------------------------------------------
-- Server version	10.6.12-MariaDB-0ubuntu0.22.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `demultiplex`
--

DROP TABLE IF EXISTS `demultiplex`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `demultiplex` (
  `trial` varchar(100) NOT NULL,
  `subject` varchar(100) NOT NULL,
  `sample` varchar(100) NOT NULL,
  `replicate` tinyint(3) unsigned NOT NULL,
  `refGenome` varchar(50) NOT NULL,
  `vectorFastaFile` varchar(100) DEFAULT NULL,
  `flags` varchar(50) DEFAULT NULL,
  `index1Seq` varchar(100) DEFAULT NULL,
  `adriftReadLinkerSeq` varchar(100) DEFAULT NULL,
  `leaderSeqHMM` varchar(100) DEFAULT NULL,
  `seqRunID` varchar(100) NOT NULL,
  PRIMARY KEY (`trial`,`subject`,`sample`,`replicate`,`refGenome`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `fragments`
--

DROP TABLE IF EXISTS `fragments`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `fragments` (
  `trial` varchar(100) NOT NULL,
  `subject` varchar(100) NOT NULL,
  `sample` varchar(100) NOT NULL,
  `replicate` tinyint(3) unsigned NOT NULL,
  `refGenome` varchar(50) NOT NULL,
  `vectorFastaFile` varchar(100) DEFAULT NULL,
  `flags` varchar(50) DEFAULT NULL,
  `seqRunID` varchar(50) DEFAULT NULL,
  `nDuplicateReads` int(10) unsigned DEFAULT NULL,
  `readID` varchar(100) DEFAULT NULL,
  `chromosome` varchar(25) DEFAULT NULL,
  `strand` varchar(1) DEFAULT NULL,
  `fragStart` int(10) unsigned DEFAULT NULL,
  `fragEnd` int(10) unsigned DEFAULT NULL,
  `leaderSeq` varchar(1000) DEFAULT NULL,
  `randomLinkerSeq` varchar(50) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `multihits`
--

DROP TABLE IF EXISTS `multihits`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `multihits` (
  `trial` varchar(100) NOT NULL,
  `subject` varchar(100) NOT NULL,
  `sample` varchar(100) NOT NULL,
  `refGenome` varchar(50) NOT NULL,
  `data` longblob DEFAULT NULL,
  PRIMARY KEY (`trial`,`subject`,`sample`,`refGenome`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sites`
--

DROP TABLE IF EXISTS `sites`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sites` (
  `trial` varchar(100) NOT NULL,
  `subject` varchar(100) NOT NULL,
  `sample` varchar(100) NOT NULL,
  `refGenome` varchar(50) NOT NULL,
  `vectorFastaFile` varchar(100) DEFAULT NULL,
  `flags` varchar(50) DEFAULT NULL,
  `posid` varchar(25) DEFAULT NULL,
  `UMIs` int(10) unsigned DEFAULT NULL,
  `sonicLengths` int(10) unsigned DEFAULT NULL,
  `reads` int(10) unsigned DEFAULT NULL,
  `repLeaderSeq` varchar(500) DEFAULT NULL,
  `replicate` tinyint(3) unsigned DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2023-10-26 17:29:33
