
echo "Diff from project_template..."
diff ../../project_templates/bwa/bwa-mem.xml  bwa-mem_v1.xml

echo "Diff from 1 to 2 (adding test)"
diff bwa-mem_v1.xml bwa-mem_v2.xml

echo "Diff from 2 to 3 (filling in test)"
diff bwa-mem_v2.xml bwa-mem_v3.xml

echo "Diff from 3 to 4 (scoring test)"
diff bwa-mem_v3.xml bwa-mem_v4.xml

echo "Diff from 4 to 5 (filling in scoring test)"
diff bwa-mem_v4.xml bwa-mem_v5.xml
