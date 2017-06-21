
echo "seqtk v1 to v2"
diff seqtk_seq_v1.xml seqtk_seq_v2.xml

echo "seqtk v2 to v3"
diff seqtk_seq_v2.xml seqtk_seq_v3.xml

echo "seqtk v3 to v4"
diff seqtk_seq_v3.xml seqtk_seq_v4.xml

echo "seqtk v4 to v5"
diff seqtk_seq_v4.xml seqtk_seq_v5.xml

echo "seqtk v5 to v6"
diff seqtk_seq_v5.xml seqtk_seq_v6.xml


echo "Diff bwa from project_template..."
diff ../../project_templates/bwa/bwa-mem.xml  bwa-mem_v1.xml

echo "Diff from 1 to 2 (adding test)"
diff bwa-mem_v1.xml bwa-mem_v2.xml

echo "Diff from 2 to 3 (filling in test)"
diff bwa-mem_v2.xml bwa-mem_v3.xml

echo "Diff from 3 to 4 (scoring test)"
diff bwa-mem_v3.xml bwa-mem_v4.xml

echo "Diff from 4 to 5 (filling in scoring test)"
diff bwa-mem_v4.xml bwa-mem_v5.xml

echo "pear from question to answer"
diff ../../project_templates/conda_exercises/exercise_1/pear.xml ../../project_templates/conda_answers/exercise_1/pear.xml

echo "fleeqtk from question to answer"
diff ../../project_templates/conda_exercises/exercise_2/fleeqtk_seq.xml ../../project_templates/conda_answers/exercise_2/fleeqtk_seq.xml
