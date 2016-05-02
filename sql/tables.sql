DROP TABLE hepBounds;
CREATE TABLE hepBounds (
  sqrtS DOUBLE NOT NULL ,
  mi INT NOT NULL CHECK (mi >=0 AND mi < 10),
  x0 DOUBLE,
  x1 DOUBLE,
  PRIMARY KEY (sqrtS, mi)
);

DROP TABLE histogramms;
CREATE TABLE histogramms (
  sqrtS DOUBLE NOT NULL ,
  mi INT NOT NULL ,
  hCol INT NOT NULL ,
  x DOUBLE,
  points INT,
  coeff DOUBLE,
  PRIMARY KEY (sqrtS, mi, hCol)
);

INSERT INTO histogramms VALUES(100,1,3,0.5,0,2);
SELECT * FROM histogramms;
DELETE FROM histogramms WHERE 1=1;

DELETE FROM hepBounds WHERE 1=1;


SELECT * FROM hepBounds;
SELECT x0,x1 FROM hepBounds WHERE mi = 1;

SELECT coeff FROM histogramms WHERE sqrtS = 100 AND mi = 0 ORDER BY hCol;
SELECT max(points) FROM histogramms WHERE sqrtS = 100 AND mi = 2;